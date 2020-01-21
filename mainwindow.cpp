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


void MainWindow::on_pushButton_compute_clicked()
{
    bool ok;

    initial_values[0] = ui->lineEdit_nu_1_0->text().toDouble(&ok);
    if (!ok)
        return;

    initial_values[1] = ui->lineEdit_nu_2_0->text().toDouble(&ok);
    if (!ok)
        return;

    initial_values[2] = ui->lineEdit_nu_3_0->text().toDouble(&ok);
    if (!ok)
        return;

    initial_values[3] = 0;

    initial_values[4] = 0;

    initial_values[5] = 0;

    final_values[0] = ui->lineEdit_nu_1_T->text().toDouble(&ok);
    if (!ok)
        return;

    final_values[1] = ui->lineEdit_nu_2_T->text().toDouble(&ok);
    if (!ok)
        return;

    final_values[2] = ui->lineEdit_nu_3_T->text().toDouble(&ok);
    if (!ok)
        return;

    final_values[3] = ui->lineEdit_x_T->text().toDouble(&ok);
    if (!ok)
        return;

    final_values[4] = ui->lineEdit_y_T->text().toDouble(&ok);
    if (!ok)
        return;

    final_values[5] = ui->lineEdit_theta_T->text().toDouble(&ok);
    if (!ok)
        return;

    t_sw = ui->lineEdit_t_sw->text().toDouble(&ok);
    if (!ok || t_sw <= 0)
        return;

    T = ui->lineEdit_T->text().toDouble(&ok);
    if (!ok || T <= t_sw)
        return;

    if (plotted)
    {
        t.clear();
        nu_1.clear();
        nu_2.clear();
        nu_3.clear();
        x.clear();
        y.clear();
        theta.clear();

        t_err.clear();
        nu_1_err.clear();
        nu_2_err.clear();
        nu_3_err.clear();
        x_err.clear();
        y_err.clear();
        theta_err.clear();
    }

    Vector<6> control = predict_control(t_sw, T, initial_values[0], final_values[0],
            initial_values[1], final_values[1], initial_values[2], final_values[2],
            final_values[3], final_values[4], final_values[5]);

    control = custom_control_find(control, t_sw, T, initial_values, final_values);

    DOPRI8_final_plot (0, T, initial_values, {control[0], control[1], control[2]},
                        {control[3], control[4], control[5]},
                        t_sw, t, nu_1, nu_2, nu_3, x, y, theta);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, 1};
    std::function<Vector<3>()> error_function = [&]() {
        static int count = 0;
        count += 1;

        if (count % 100 == 0)
            return Vector<3>({d(gen), d(gen), d(gen)});
        else
            return Vector<3>({0, 0, 0});
    };

    DOPRI8_final_plot (0, T, initial_values, {control[0], control[1], control[2]},
                        {control[3], control[4], control[5]},
                        t_sw, t_err, nu_1_err, nu_2_err, nu_3_err, x_err, y_err, theta_err, error_function);

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
    }


    ui->textBrowser_controls->setText("U_expected = " + QString::number(control[0], 'g', 6)
            + " " + QString::number(control[1], 'g', 6) + " " + QString::number(control[2], 'g', 6)
            + " " + QString::number(control[3], 'g', 6) + " " + QString::number(control[4], 'g', 6)
            + " " + QString::number(control[5], 'g', 6));

    trajectory_minus = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_minus_err = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus_err = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_minus, data_plus, data_minus_err, data_plus_err;

    QPen pen_minus(Qt::gray);
    QPen pen_plus(Qt::yellow);
    trajectory_minus->setPen(pen_minus);
    trajectory_plus->setPen(pen_plus);

    QPen pen_minus_err(Qt::DashLine);
    pen_minus_err.setColor(Qt::darkRed);
    QPen pen_plus_err(Qt::DashLine);
    pen_minus_err.setColor(Qt::darkCyan);
    trajectory_minus_err->setPen(pen_minus_err);
    trajectory_plus_err->setPen(pen_plus_err);

    int i = 0;

    for (i = 0; t[i] < t_sw; i++)
        data_minus.append(QCPCurveData(i, x[i], y[i]));

    for (; i < x.length(); i++)
        data_plus.append(QCPCurveData(i, x[i], y[i]));

    trajectory_minus->data()->set(data_minus, true);
    trajectory_plus->data()->set(data_plus, true);

    for (i = 0; t_err[i] < t_sw; i++)
        data_minus_err.append(QCPCurveData(i, x_err[i], y_err[i]));

    for (; i < x_err.length(); i++)
        data_plus_err.append(QCPCurveData(i, x_err[i], y_err[i]));

    trajectory_minus_err->data()->set(data_minus_err, true);
    trajectory_plus_err->data()->set(data_plus_err, true);

    double x_max_1 = *std::max_element(x.begin(), x.end());
    double x_min_1 = *std::min_element(x.begin(), x.end());
    double y_max_1 = *std::max_element(y.begin(), y.end());
    double y_min_1 = *std::min_element(y.begin(), y.end());

    double x_max_2 = *std::max_element(x_err.begin(), x_err.end());
    double x_min_2 = *std::min_element(x_err.begin(), x_err.end());
    double y_max_2 = *std::max_element(y_err.begin(), y_err.end());
    double y_min_2 = *std::min_element(y_err.begin(), y_err.end());

    double x_max = std::max(x_max_1, x_max_2);
    double x_min = std::min(x_min_1, x_min_2);
    double y_max = std::max(y_max_1, y_max_2);
    double y_min = std::min(y_min_1, y_min_2);

    ui->PlotWidget_trajectory->xAxis->setRange(x_min - (x_max - x_min) * 0.05, x_max + (x_max - x_min) * 0.05);
    ui->PlotWidget_trajectory->yAxis->setRange(y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05);
    ui->PlotWidget_trajectory->xAxis->setLabel("x");
    ui->PlotWidget_trajectory->yAxis->setLabel("y");

    ui->PlotWidget_trajectory->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_trajectory->replot();

    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/trajectory_t_sw_"
                                + QString::number(t_sw, 'g', 4) + "_T_" + QString::number(T, 'g', 4)
                                + "_nu_1_0_" + QString::number(initial_values[0], 'g', 4)
                                + "_nu_2_0_" + QString::number(initial_values[1], 'g', 4)
                                + "_nu_3_0_" + QString::number(initial_values[2], 'g', 4)
                                + "_nu_1_T_" + QString::number(final_values[0], 'g', 4)
                                + "_nu_2_T_" + QString::number(final_values[1], 'g', 4)
                                + "_nu_3_T_" + QString::number(final_values[2], 'g', 4)
                                + "_x_T_" + QString::number(final_values[3], 'g', 4)
                                + "_y_T_" + QString::number(final_values[4], 'g', 4)
                                + "_theta_T_" + QString::number(final_values[5], 'g', 4)
                                + "_control_error.pdf");

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
