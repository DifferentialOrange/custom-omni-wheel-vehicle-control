#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    plotted(false)
{
    ui->setupUi(this);

    ui->lineEdit_nu_1_0->setText("0");
    ui->lineEdit_nu_2_0->setText("0");
    ui->lineEdit_nu_3_0->setText("0");
    ui->lineEdit_nu_1_T->setText("0");
    ui->lineEdit_nu_2_T->setText("0");
    ui->lineEdit_nu_3_T->setText("0");
    ui->lineEdit_x_T->setText("10");
    ui->lineEdit_y_T->setText("10");
    ui->lineEdit_theta_T->setText("0");
    ui->lineEdit_t_sw->setText("5");
    ui->lineEdit_T->setText("10");
}

MainWindow::~MainWindow()
{
    delete ui;
}

static double energy(Vector<6> control, double t_sw, double T)
{
    return t_sw * (control[0] * control[0] + control[1] * control[1] + control[2] * control[2]) + \
            (T - t_sw) * (control[3] * control[3] + control[4] * control[4] + control[5] * control[5]);
}

void MainWindow::on_pushButton_compute_clicked()
{
    bool ok;

    initial_values[0] = 0;

    initial_values[1] = 0;

    initial_values[2] = 0;

    initial_values[3] = 0;

    initial_values[4] = 0;

    initial_values[5] = 0;

    final_values[0] = 0;

    final_values[1] = 0;

    final_values[2] = 0;

    final_values[3] = ui->lineEdit_x_T->text().toDouble(&ok);
    if (!ok)
        return;

    final_values[4] = ui->lineEdit_y_T->text().toDouble(&ok);
    if (!ok)
        return;

    final_values[5] = ui->lineEdit_theta_T->text().toDouble(&ok);
    if (!ok)
        return;

    T = ui->lineEdit_T->text().toDouble(&ok);
    if (!ok)
        return;

    double min_energy_dynamics = 1e30;

    for (double t_sw = T / 1000; t_sw <= T * 999 / 1000; t_sw += T / 1000)
    {
        Vector<6> control = predict_control(t_sw, T, initial_values[0], final_values[0],
                initial_values[1], final_values[1], initial_values[2], final_values[2],
                final_values[3], final_values[4], final_values[5]);

        qDebug() << "Computed dynamic strategy energy for t_sw = " << t_sw << "\n";
        min_energy_dynamics = std::min(min_energy_dynamics, energy(control, t_sw, T));
    }

    QVector<double> s_vect;
    QVector<double> energy_trajectory;
    QVector<double> energy_dynamics;

    for (double s = T / 100; s <= T * 99 / 100; s += T / 100)
    {
        double min_energy_trajectory_turn = 1e30;
        for (double t_sw_turn = s / 10; t_sw_turn < s * 9 / 10; t_sw_turn += T / 100)
        {
            Vector<6> control = predict_control(t_sw_turn, s, initial_values[0], final_values[0],
                    initial_values[1], final_values[1], initial_values[2], final_values[2],
                    0, 0, final_values[5]);

            qDebug() << "Computed turn energy for t_sw = " << t_sw_turn << "\n";
            min_energy_trajectory_turn = std::min(min_energy_trajectory_turn, energy(control, t_sw_turn, s));
        }

        double min_energy_trajectory_line = 1e30;
        for (double t_sw_line = (T - s) / 10; t_sw_line < (T - s) * 9 / 10; t_sw_line += T / 100)
        {
            Vector<6> control = predict_control(t_sw_line, T - s, initial_values[0], final_values[0],
                    initial_values[1], final_values[1], initial_values[2], final_values[2],
                    final_values[3], final_values[4], 0);

            qDebug() << "Computed line movement energy for t_sw = " << t_sw_line << "\n";
            min_energy_trajectory_line = std::min(min_energy_trajectory_line, energy(control, t_sw_line, T - s));
        }

        qDebug() << "Computed energy for s = " << s << " is" << min_energy_trajectory_turn + min_energy_trajectory_line << "\n";
        s_vect.push_back(s);
        energy_trajectory.push_back(min_energy_trajectory_turn + min_energy_trajectory_line);
        energy_dynamics.push_back(min_energy_dynamics);
    }

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
    }

    ui->PlotWidget_trajectory->addGraph();
    ui->PlotWidget_trajectory->graph(0)->setPen(QPen(Qt::green));
    ui->PlotWidget_trajectory->graph(0)->setData(s_vect, energy_trajectory);
    ui->PlotWidget_trajectory->addGraph();
    ui->PlotWidget_trajectory->graph(1)->setPen(QPen(Qt::red));
    ui->PlotWidget_trajectory->graph(1)->setData(s_vect, energy_dynamics);

    ui->PlotWidget_trajectory->graph(0)->rescaleAxes();
    ui->PlotWidget_trajectory->graph(1)->rescaleAxes(true);

    ui->PlotWidget_trajectory->xAxis->setLabel("s");
    ui->PlotWidget_trajectory->yAxis->setLabel("E");

    ui->PlotWidget_trajectory->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_trajectory->replot();

    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/energy_comparison_T_"
                                + QString::number(T, 'g', 4)
                                + "_x_T_" + QString::number(final_values[3], 'g', 4)
                                + "_y_T_" + QString::number(final_values[4], 'g', 4)
                                + "_theta_T_" + QString::number(final_values[5], 'g', 4)
                                + ".pdf");

    plotted = true;
}

void MainWindow::on_pushButton_generate_clicked()
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> nu_1_2(1.0, -1.0);
    std::uniform_real_distribution<> nu_3(0.2, -0.2);
    std::uniform_real_distribution<> x_y(-50, 50);
    std::uniform_real_distribution<> theta(- 4 * M_PI, 4 * M_PI);
    std::uniform_real_distribution<> time(1, 50);

    double t_1 = time(gen), t_2 = time(gen);

    if (t_1 > t_2)
    {
        ui->lineEdit_T->setText(QString::number(t_1));
        ui->lineEdit_t_sw->setText(QString::number(t_2));
    }
    else
    {
        ui->lineEdit_T->setText(QString::number(t_2));
        ui->lineEdit_t_sw->setText(QString::number(t_1));
    }

    ui->lineEdit_nu_1_0->setText(QString::number(nu_1_2(gen)));
    ui->lineEdit_nu_2_0->setText(QString::number(nu_1_2(gen)));
    ui->lineEdit_nu_1_T->setText(QString::number(nu_1_2(gen)));
    ui->lineEdit_nu_2_T->setText(QString::number(nu_1_2(gen)));

    ui->lineEdit_nu_3_0->setText(QString::number(nu_3(gen)));
    ui->lineEdit_nu_3_T->setText(QString::number(nu_3(gen)));

    ui->lineEdit_x_T->setText(QString::number(x_y(gen)));
    ui->lineEdit_y_T->setText(QString::number(x_y(gen)));
    ui->lineEdit_theta_T->setText(QString::number(theta(gen)));
}
