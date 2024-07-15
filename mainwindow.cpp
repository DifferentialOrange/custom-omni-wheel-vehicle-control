#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include.h"
#include <float.h>

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

double e_c(double t_sw, double T) {
    double c = parameters::symmetrical::kappa;
    return ((exp(c * t_sw) - exp(c * T)) * (exp(c * t_sw) - exp(c * T)) * t_sw + \
            (exp(c * t_sw) - 1) * (exp(c * t_sw) - 1) * (T - t_sw)) /
            (t_sw - t_sw * exp(c * T) - T + exp(c * t_sw)) /
            (t_sw - t_sw * exp(c * T) - T + exp(c * t_sw));
}

double e_c_3(double t_sw, double T) {
    double c_3 = parameters::symmetrical::kappa3;
    return ((exp(c_3 * t_sw) - exp(c_3 * T)) * (exp(c_3 * t_sw) - exp(c_3 * T)) * t_sw + \
            (exp(c_3 * t_sw) - 1) * (exp(c_3 * t_sw) - 1) * (T - t_sw)) /
            (t_sw - t_sw * exp(c_3 * T) - T + exp(c_3 * t_sw)) /
            (t_sw - t_sw * exp(c_3 * T) - T + exp(c_3 * t_sw));
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
        t_symm.clear();
        nu_1_symm.clear();
        nu_2_symm.clear();
        nu_3_symm.clear();
        x_symm.clear();
        y_symm.clear();
        theta_symm.clear();
    }


    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
    }

    QVector<double> tv, detv;

    T = 1;

    int i = 0;
//    tv.append(0);
//    detv.append(0);
    for (double t_sw = 0.01; t_sw < T - 0.01; t_sw += 0.001, ++i) {
        double det = e_c(t_sw, T);
        tv.append(t_sw);
        detv.append(det);
    }
//    tv.append(T);
//    detv.append(0);

    double minv = DBL_MAX;
    double minarg = -1;
    for (double t_sw = 0.49; t_sw < 0.51; t_sw += 1e-9) {
        double v = e_c(t_sw, T);
        if (v < minv) {
            minarg = t_sw;
            minv = v;
        }
    }
    qDebug() << "min" << minarg << '\n';

    ui->PlotWidget_trajectory->addGraph();
    ui->PlotWidget_trajectory->graph(0)->setData(tv, detv);

    QPen pen(Qt::magenta);
    pen.setWidth(3);
    ui->PlotWidget_trajectory->graph(0)->setPen(pen);

    ui->PlotWidget_trajectory->xAxis->setLabel("t_sw");
    ui->PlotWidget_trajectory->yAxis->setLabel("e_c");
    ui->PlotWidget_trajectory->rescaleAxes();
    ui->PlotWidget_trajectory->xAxis->setRange(0, T);

    ui->PlotWidget_trajectory->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_trajectory->replot();

    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/e_c_T_"
                                + QString::number(T, 'g', 4)
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
