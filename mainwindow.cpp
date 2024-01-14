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
//    ui->lineEdit_t_1->setText("5");
//    ui->lineEdit_t_2->setText("10");
//    ui->lineEdit_U_1_1->setText("80");
//    ui->lineEdit_U_2_1->setText("80");
//    ui->lineEdit_U_3_1->setText("-70");
//    ui->lineEdit_U_1_2->setText("0");
//    ui->lineEdit_U_2_2->setText("0");
//    ui->lineEdit_U_3_2->setText("0");
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

    Vector<3> control_1, control_2;

    control_1[0] = ui->lineEdit_U_1_m->text().toDouble(&ok);
    if (!ok)
        return;

    control_1[1] = ui->lineEdit_U_2_m->text().toDouble(&ok);
    if (!ok)
        return;

    control_1[2] = ui->lineEdit_U_3_m->text().toDouble(&ok);
    if (!ok)
        return;

    control_2[0] = ui->lineEdit_U_1_p->text().toDouble(&ok);
    if (!ok)
        return;

    control_2[1] = ui->lineEdit_U_2_p->text().toDouble(&ok);
    if (!ok)
        return;

    control_2[2] = ui->lineEdit_U_3_p->text().toDouble(&ok);
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

        P_real.clear();
        P_advice.clear();

        N_1.clear();
        N_2.clear();
        N_3.clear();

        t.clear();
        nu_1.clear();
        nu_2.clear();
        nu_3.clear();
        x.clear();
        y.clear();
        theta.clear();
    }

    DOPRI8_symmetrical_plot (0, T, initial_values, control_1,
                             control_2, t_sw,
                             t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                             x_symm, y_symm, theta_symm,
                             P_real, P_advice, N_1, N_2, N_3);

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_P->clearPlottables();
        ui->PlotWidget_N->clearPlottables();
    }

    trajectory_minus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_minus_symm, data_plus_symm;

    QPen pen_minus_symm(Qt::blue);
    QPen pen_plus_symm(Qt::magenta);
    trajectory_minus_symm->setPen(pen_minus_symm);
    trajectory_plus_symm->setPen(pen_plus_symm);

    int i = 0;
    int i_boundary;
    for (i = 0; t_symm[i] < t_sw; i++)
        data_minus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

    i_boundary = i;

    for (; i < x_symm.length(); i++)
        data_plus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

    trajectory_minus_symm->data()->set(data_minus_symm, true);
    trajectory_plus_symm->data()->set(data_plus_symm, true);

    double g = 9.81;
    double bad_advice_1 = (-3 * sqrt(3) * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * \
                           nu_1_symm[i_boundary] * nu_3_symm[i_boundary] - \
                           3 * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * nu_2_symm[i_boundary] * nu_3_symm[i_boundary] - \
                            3 * parameters::c2 * nu_1_symm[i_boundary] / 2 + 3 * sqrt(3) * parameters::c2 * nu_2_symm[i_boundary] /2 - \
                            g*parameters::symmetrical::rho) / parameters::c1;

    double bad_advice_2 = (3 * parameters::lambda * parameters::lambda / parameters::symmetrical::Lambda * nu_2_symm[i_boundary] * nu_3_symm[i_boundary] + \
                            3 * parameters::c2 * nu_1_symm[i_boundary] - \
                            g*parameters::symmetrical::rho) / parameters::c1;

    double bad_advice_3 = (3 * sqrt(3) * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * \
                           nu_1_symm[i_boundary] * nu_3_symm[i_boundary] - \
                           3 * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * nu_2_symm[i_boundary] * nu_3_symm[i_boundary] - \
                            3 * parameters::c2 * nu_1_symm[0] / 2 - 3 * sqrt(3) * parameters::c2 * nu_2_symm[i_boundary] /2 - \
                            g*parameters::symmetrical::rho) / parameters::c1;

    qDebug() << "try U_1_plus + U_2_plus - 2 U_3_plus > " << bad_advice_1 << " for more fun result" << '\n';
    qDebug() << "try U_1_plus - 2 U_2_plus + U_3_plus > " << bad_advice_2 << " for more fun result" << '\n';
    qDebug() << "try -2 U_1_plus + U_2_plus + U_3_plus > " << bad_advice_3 << " for more fun result" << '\n';

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
                                + ".pdf");
    qDebug() << "PlotWidget_trajectory" << "\n";

    double P_max_1 = *std::max_element(P_real.begin(), P_real.end());
    double P_max_2 = *std::max_element(P_advice.begin(), P_advice.end());
    double P_max = std::max(P_max_1, P_max_2);

    double P_min_1 = *std::min_element(P_real.begin(), P_real.end());
    double P_min_2 = *std::min_element(P_advice.begin(), P_advice.end());
    double P_min = std::min(P_min_1, P_min_2);

    ui->PlotWidget_P->legend->setVisible(true);

    QPen pen_advice(Qt::green);
    pen_advice.setStyle(Qt::DashLine);

    ui->PlotWidget_P->addGraph();
    ui->PlotWidget_P->graph(0)->setData(t_symm, P_advice);
    ui->PlotWidget_P->graph(0)->setName("advice");
    ui->PlotWidget_P->graph(0)->setPen(pen_advice);

    ui->PlotWidget_P->addGraph();
    ui->PlotWidget_P->graph(1)->setData(t_symm, P_real);
    ui->PlotWidget_P->graph(1)->setName("real");
    ui->PlotWidget_P->graph(1)->setPen(QPen(Qt::red));

    ui->PlotWidget_P->xAxis->setRange(0, t_symm.last());
    ui->PlotWidget_P->yAxis->setRange(- 0.05 + P_min, P_max + 0.05);
    ui->PlotWidget_P->xAxis->setLabel("t");
    ui->PlotWidget_P->yAxis->setLabel("P");
    ui->PlotWidget_P->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_P->replot();

    ui->PlotWidget_P->savePdf("../custom-omni-wheel-vehicle-control/PICS/power_"
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
                                + ".pdf");

    double N_max_1 = *std::max_element(N_1.begin(), N_1.end());
    double N_max_2 = *std::max_element(N_2.begin(), N_2.end());
    double N_max_3 = *std::max_element(N_3.begin(), N_3.end());
    double N_min_1 = *std::min_element(N_1.begin(), N_1.end());
    double N_min_2 = *std::min_element(N_2.begin(), N_2.end());
    double N_min_3 = *std::min_element(N_3.begin(), N_3.end());
    double N_max = std::max(std::max(N_max_1, N_max_2), N_max_3);
    double N_min = std::min(std::min(N_min_1, N_min_2), N_min_3);

    ui->PlotWidget_N->legend->setVisible(true);

    ui->PlotWidget_N->addGraph();
    ui->PlotWidget_N->graph(0)->setData(t_symm, N_1);
    ui->PlotWidget_N->graph(0)->setName("N_1");
    ui->PlotWidget_N->graph(0)->setPen(QPen(Qt::blue));

    ui->PlotWidget_N->addGraph();
    ui->PlotWidget_N->graph(1)->setData(t_symm, N_2);
    ui->PlotWidget_N->graph(1)->setName("N_2");
    ui->PlotWidget_N->graph(1)->setPen(QPen(Qt::magenta));

    ui->PlotWidget_N->addGraph();
    ui->PlotWidget_N->graph(2)->setData(t_symm, N_3);
    ui->PlotWidget_N->graph(2)->setName("N_3");
    ui->PlotWidget_N->graph(2)->setPen(QPen(Qt::cyan));

    ui->PlotWidget_N->xAxis->setRange(0, t_symm.last());
    ui->PlotWidget_N->yAxis->setRange(N_min - 0.05, N_max + 0.05);
    ui->PlotWidget_N->xAxis->setLabel("t");
    ui->PlotWidget_N->yAxis->setLabel("N");
    ui->PlotWidget_N->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_N->replot();

    ui->PlotWidget_N->savePdf("../custom-omni-wheel-vehicle-control/PICS/N_"
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

//    ui->lineEdit_nu_1_0->setText(QString::number(nu_1_2(gen)));
//    ui->lineEdit_nu_2_0->setText(QString::number(nu_1_2(gen)));
//    ui->lineEdit_nu_1_T->setText(QString::number(nu_1_2(gen)));
//    ui->lineEdit_nu_2_T->setText(QString::number(nu_1_2(gen)));

//    ui->lineEdit_nu_3_0->setText(QString::number(nu_3(gen)));
//    ui->lineEdit_nu_3_T->setText(QString::number(nu_3(gen)));

//    ui->lineEdit_x_T->setText(QString::number(x_y(gen)));
//    ui->lineEdit_y_T->setText(QString::number(x_y(gen)));
//    ui->lineEdit_theta_T->setText(QString::number(theta(gen)));
}
