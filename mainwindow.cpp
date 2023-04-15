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
    ui->lineEdit_t_1->setText("5");
    ui->lineEdit_t_2->setText("10");
    ui->lineEdit_U_1_1->setText("80");
    ui->lineEdit_U_2_1->setText("80");
    ui->lineEdit_U_3_1->setText("-70");
    ui->lineEdit_U_1_2->setText("0");
    ui->lineEdit_U_2_2->setText("0");
    ui->lineEdit_U_3_2->setText("0");
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_compute_clicked()
{
    bool ok;

    initial_values_1[0] = ui->lineEdit_nu_1_0->text().toDouble(&ok);
    if (!ok)
        return;

    initial_values_1[1] = ui->lineEdit_nu_2_0->text().toDouble(&ok);
    if (!ok)
        return;

    initial_values_1[2] = ui->lineEdit_nu_3_0->text().toDouble(&ok);
    if (!ok)
        return;

    initial_values_1[3] = 0;

    initial_values_1[4] = 0;

    initial_values_1[5] = 0;

    Vector<3> control_1, control_2;

    control_1[0] = ui->lineEdit_U_1_1->text().toDouble(&ok);
    if (!ok)
        return;

    control_1[1] = ui->lineEdit_U_2_1->text().toDouble(&ok);
    if (!ok)
        return;

    control_1[2] = ui->lineEdit_U_3_1->text().toDouble(&ok);
    if (!ok)
        return;

    control_2[0] = ui->lineEdit_U_1_2->text().toDouble(&ok);
    if (!ok)
        return;

    control_2[1] = ui->lineEdit_U_2_2->text().toDouble(&ok);
    if (!ok)
        return;

    control_2[2] = ui->lineEdit_U_3_2->text().toDouble(&ok);
    if (!ok)
        return;

    double t_1 = ui->lineEdit_t_1->text().toDouble(&ok);
    if (!ok || t_1 <= 0)
        return;

    double t_2 = ui->lineEdit_t_2->text().toDouble(&ok);
    if (!ok || t_2 <= t_1)
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

        N_symm_1.clear();
        N_symm_2.clear();
        N_symm_3.clear();

    }

    DOPRI8_symmetrical_plot (0, t_1, initial_values_1, control_1,
                             t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                             x_symm, y_symm, theta_symm,
                             N_symm_1, N_symm_2, N_symm_3);

    qDebug() << "plot 1 ok" << '\n';

    initial_values_2[0] = nu_1_symm.takeLast();

    initial_values_2[1] = nu_2_symm.takeLast();

    initial_values_2[2] = nu_3_symm.takeLast();

    initial_values_2[3] = x_symm.takeLast();

    initial_values_2[4] = y_symm.takeLast();

    initial_values_2[5] = theta_symm.takeLast();

    double g = 9.81;
    double bad_advice_1 = (-3 * sqrt(3) * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * \
                           initial_values_2[0] * initial_values_2[2] - \
                           3 * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * initial_values_2[1] * initial_values_2[2] - \
                            3 * parameters::c2 * initial_values_2[0] / 2 + 3 * sqrt(3) * parameters::c2 * initial_values_2[1] /2 - \
                            g*parameters::symmetrical::rho) / parameters::c1;

    double bad_advice_2 = (3 * parameters::lambda * parameters::lambda / parameters::symmetrical::Lambda * initial_values_2[1] * initial_values_2[2] + \
                            3 * parameters::c2 * initial_values_2[0] - \
                            g*parameters::symmetrical::rho) / parameters::c1;

    double bad_advice_3 = (3 * sqrt(3) * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * \
                           initial_values_2[0] * initial_values_2[2] - \
                           3 * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * initial_values_2[1] * initial_values_2[2] - \
                            3 * parameters::c2 * initial_values_2[0] / 2 - 3 * sqrt(3) * parameters::c2 * initial_values_2[1] /2 - \
                            g*parameters::symmetrical::rho) / parameters::c1;

    qDebug() << "try    U_1 +   U_2 - 2 U_3 > " << bad_advice_1 << " for more fun result" << '\n';
    qDebug() << "try    U_1 - 2 U_2 +   U_3 > " << bad_advice_2 << " for more fun result" << '\n';
    qDebug() << "try -2 U_1 +   U_2 +   U_3 > " << bad_advice_3 << " for more fun result" << '\n';

    double h = (1 + 3.0 * parameters::lambda * parameters::lambda / 2) * (initial_values_2[0] * initial_values_2[0] + initial_values_2[1] * initial_values_2[1]) +
            (1 + 3.0 * parameters::symmetrical::rho * parameters::symmetrical::rho * parameters::lambda * parameters::lambda /
             parameters::symmetrical::Lambda / parameters::symmetrical::Lambda) * initial_values_2[2] * initial_values_2[2];

    double good_advice = g * parameters::symmetrical::rho / parameters::c1 / sqrt(6) -
            (3 + sqrt(3)) * (parameters::lambda * parameters::lambda * h / sqrt(parameters::symmetrical::Lambda * parameters::symmetrical::Lambda +
                                                                                3 * parameters::symmetrical::rho * parameters::symmetrical::rho *
                                                                                parameters::lambda * parameters::lambda) +
                             parameters::c2 * sqrt(h)) / 2 / parameters::c1 / sqrt(2 + 3 * parameters::lambda * parameters::lambda);

    qDebug() << "keep U_1^2 + U_2^2 + U_3^2 < " << good_advice << " for satisfying result" << '\n';


    DOPRI8_symmetrical_plot (t_1, t_2, initial_values_2, control_2,
                             t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                             x_symm, y_symm, theta_symm,
                             N_symm_1, N_symm_2, N_symm_3);

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_N_symm->clearPlottables();
    }

//    TO DO: show error in text panel

    trajectory_minus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_minus_symm, data_plus_symm;

    QPen pen_minus_symm(Qt::blue);
    QPen pen_plus_symm(Qt::magenta);
    trajectory_minus_symm->setPen(pen_minus_symm);
    trajectory_plus_symm->setPen(pen_plus_symm);

    int i = 0;
    for (i = 0; t_symm[i] < t_1; i++)
        data_minus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

    for (; i < x_symm.length(); i++)
        data_plus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

    trajectory_minus_symm->data()->set(data_minus_symm, true);
    trajectory_plus_symm->data()->set(data_plus_symm, true);


    double x_max = *std::max_element(x_symm.begin(), x_symm.end());
    double x_min = *std::min_element(x_symm.begin(), x_symm.end());
    double y_max = *std::max_element(y_symm.begin(), y_symm.end());
    double y_min = *std::min_element(y_symm.begin(), y_symm.end());

    ui->PlotWidget_trajectory->xAxis->setRange(x_min - (x_max - x_min) * 0.05, x_max + (x_max - x_min) * 0.05);
    ui->PlotWidget_trajectory->yAxis->setRange(y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05);
    ui->PlotWidget_trajectory->xAxis->setLabel("x");
    ui->PlotWidget_trajectory->yAxis->setLabel("y");

    ui->PlotWidget_trajectory->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_trajectory->replot();

//    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/trajectory_friction_t_sw_"
//                                + QString::number(t_sw, 'g', 4) + "_T_" + QString::number(T, 'g', 4)
//                                + "_nu_1_0_" + QString::number(initial_values[0], 'g', 4)
//                                + "_nu_2_0_" + QString::number(initial_values[1], 'g', 4)
//                                + "_nu_3_0_" + QString::number(initial_values[2], 'g', 4)
//                                + "_nu_1_T_" + QString::number(final_values[0], 'g', 4)
//                                + "_nu_2_T_" + QString::number(final_values[1], 'g', 4)
//                                + "_nu_3_T_" + QString::number(final_values[2], 'g', 4)
//                                + "_x_T_" + QString::number(final_values[3], 'g', 4)
//                                + "_y_T_" + QString::number(final_values[4], 'g', 4)
//                                + "_theta_T_" + QString::number(final_values[5], 'g', 4)
//                                + "_mu_n_" + QString::number(parameters::mu_n, 'g', 6)
//                                + "_mu_tau_" + QString::number(parameters::mu_tau, 'g', 6) + ".pdf");


    double N_symm_1_max = *std::max_element(N_symm_1.begin(), N_symm_1.end());
    double N_symm_2_max = *std::max_element(N_symm_2.begin(), N_symm_2.end());
    double N_symm_3_max = *std::max_element(N_symm_3.begin(), N_symm_3.end());
    double N_symm_max = std::max(N_symm_1_max, std::max(N_symm_2_max, N_symm_3_max));

    double N_symm_1_min = *std::min_element(N_symm_1.begin(), N_symm_1.end());
    double N_symm_2_min = *std::min_element(N_symm_2.begin(), N_symm_2.end());
    double N_symm_3_min = *std::min_element(N_symm_3.begin(), N_symm_3.end());
    double N_symm_min = std::min(N_symm_1_min, std::min(N_symm_2_min, N_symm_3_min));

    ui->PlotWidget_N_symm->legend->setVisible(true);

    ui->PlotWidget_N_symm->addGraph();
    ui->PlotWidget_N_symm->graph(0)->setData(t_symm, N_symm_1);
    ui->PlotWidget_N_symm->graph(0)->setName("N_1");
    ui->PlotWidget_N_symm->graph(0)->setPen(QPen(Qt::red));

    ui->PlotWidget_N_symm->addGraph();
    ui->PlotWidget_N_symm->graph(1)->setData(t_symm, N_symm_2);
    ui->PlotWidget_N_symm->graph(1)->setName("N_2");
    ui->PlotWidget_N_symm->graph(1)->setPen(QPen(Qt::green));

    ui->PlotWidget_N_symm->addGraph();
    ui->PlotWidget_N_symm->graph(2)->setData(t_symm, N_symm_3);
    ui->PlotWidget_N_symm->graph(2)->setName("N_3");
    ui->PlotWidget_N_symm->graph(2)->setPen(QPen(Qt::blue));

    ui->PlotWidget_N_symm->xAxis->setRange(0, t_2);
    ui->PlotWidget_N_symm->yAxis->setRange(N_symm_min - 0.05, N_symm_max + 0.05);
    ui->PlotWidget_N_symm->xAxis->setLabel("t");
    ui->PlotWidget_N_symm->yAxis->setLabel("N_nonholonomic");
    ui->PlotWidget_N_symm->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_N_symm->replot();

//    ui->PlotWidget_N_symm->savePdf("../custom-omni-wheel-vehicle-control/PICS/N_symm_t_sw_"
//                                + QString::number(t_sw, 'g', 4) + "_T_" + QString::number(T, 'g', 4)
//                                + "_nu_1_0_" + QString::number(initial_values[0], 'g', 4)
//                                + "_nu_2_0_" + QString::number(initial_values[1], 'g', 4)
//                                + "_nu_3_0_" + QString::number(initial_values[2], 'g', 4)
//                                + "_nu_1_T_" + QString::number(final_values[0], 'g', 4)
//                                + "_nu_2_T_" + QString::number(final_values[1], 'g', 4)
//                                + "_nu_3_T_" + QString::number(final_values[2], 'g', 4)
//                                + "_x_T_" + QString::number(final_values[3], 'g', 4)
//                                + "_y_T_" + QString::number(final_values[4], 'g', 4)
//                                + "_theta_T_" + QString::number(final_values[5], 'g', 4)
//                                + "_mu_n_" + QString::number(parameters::mu_n, 'g', 6)
//                                + "_mu_tau_" + QString::number(parameters::mu_tau, 'g', 6) + ".pdf");

    plotted = true;
}


