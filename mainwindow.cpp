#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    plotted(false)
{
    ui->setupUi(this);

//    ui->lineEdit_nu_1_0->setText("0");
//    ui->lineEdit_nu_2_0->setText("0");
//    ui->lineEdit_nu_3_0->setText("0");
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
        PT_advice.clear();

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

        U1.clear();
        U2.clear();
        U3.clear();
    }

    DOPRI8_symmetrical_plot (0, T, initial_values, control_1,
                             control_2, t_sw,
                             t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                             x_symm, y_symm, theta_symm,
                             P_real, P_advice, PT_advice, N_1, N_2, N_3, U1, U2, U3);

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_theta->clearPlottables();
        ui->PlotWidget_P->clearPlottables();
        ui->PlotWidget_PT->clearPlottables();
        ui->PlotWidget_N->clearPlottables();
        ui->PlotWidget_U->clearPlottables();
        ui->PlotWidget_nu_12->clearPlottables();
        ui->PlotWidget_nu_3->clearPlottables();
    }

    trajectory_minus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_minus_symm, data_plus_symm;

    QPen pen_minus_symm(Qt::magenta);
    QPen pen_plus_symm(Qt::magenta);
    trajectory_minus_symm->setPen(pen_minus_symm);
    trajectory_plus_symm->setPen(pen_plus_symm);

    int i = 0;
    int i_boundary;
    bool first_half = false;
    bool second_half = false;
    int current_lap = 0;
    for (i = 0; t_symm[i] < t_sw; i++) {
        data_minus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));
        if (y_symm[i] >= 0 && !first_half){
            current_lap++;
            first_half = true;
            second_half = false;
        }
        if (y_symm[i] < 0 && !second_half){
            second_half = true;
            first_half = false;
        }
    }


    i_boundary = i;

    for (; i < x_symm.length(); i++){
        data_plus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));
        if (y_symm[i] >= 0 && !first_half){
            current_lap++;
            first_half = true;
            second_half = false;
        }
        if (y_symm[i] < 0 && !second_half){
            second_half = true;
            first_half = false;
        }
    }

    qDebug() << "stopped on the " << current_lap << " lap" << '\n';
    qDebug() << "x " << x_symm[x_symm.length() - 1] << "y " << y_symm[x_symm.length() - 1] << '\n';


    trajectory_minus_symm->data()->set(data_minus_symm, true);
    trajectory_plus_symm->data()->set(data_plus_symm, true);

//    double g = 9.81;
//    double bad_advice_1 = (-3 * sqrt(3) * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * \
//                           nu_1_symm[i_boundary] * nu_3_symm[i_boundary] - \
//                           3 * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * nu_2_symm[i_boundary] * nu_3_symm[i_boundary] - \
//                            3 * parameters::c2 * nu_1_symm[i_boundary] / 2 + 3 * sqrt(3) * parameters::c2 * nu_2_symm[i_boundary] /2 - \
//                            g*parameters::symmetrical::rho) / parameters::c1;

//    double bad_advice_2 = (3 * parameters::lambda * parameters::lambda / parameters::symmetrical::Lambda * nu_2_symm[i_boundary] * nu_3_symm[i_boundary] + \
//                            3 * parameters::c2 * nu_1_symm[i_boundary] - \
//                            g*parameters::symmetrical::rho) / parameters::c1;

//    double bad_advice_3 = (3 * sqrt(3) * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * \
//                           nu_1_symm[i_boundary] * nu_3_symm[i_boundary] - \
//                           3 * parameters::lambda * parameters::lambda / 2 / parameters::symmetrical::Lambda * nu_2_symm[i_boundary] * nu_3_symm[i_boundary] - \
//                            3 * parameters::c2 * nu_1_symm[0] / 2 - 3 * sqrt(3) * parameters::c2 * nu_2_symm[i_boundary] /2 - \
//                            g*parameters::symmetrical::rho) / parameters::c1;

//    qDebug() << "try U_1_plus + U_2_plus - 2 U_3_plus > " << bad_advice_1 << " for more fun result" << '\n';
//    qDebug() << "try U_1_plus - 2 U_2_plus + U_3_plus > " << bad_advice_2 << " for more fun result" << '\n';
//    qDebug() << "try -2 U_1_plus + U_2_plus + U_3_plus > " << bad_advice_3 << " for more fun result" << '\n';

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

    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_trajectory_explicit.pdf");

    double P_max_1 = *std::max_element(P_real.begin(), P_real.end());
    double P_max_2 = *std::max_element(P_advice.begin(), P_advice.end());
    double P_max_3 = *std::max_element(PT_advice.begin(), PT_advice.end());
    double P_max1 = std::max(P_max_1, P_max_2);
    double P_max = std::max(P_max1, P_max_3);

    ui->PlotWidget_P->legend->setVisible(true);

    QPen pen_advice(Qt::red);
    pen_advice.setStyle(Qt::DashLine);
    ui->PlotWidget_P->addGraph();
    ui->PlotWidget_P->graph(0)->setData(t_symm, P_advice);
    ui->PlotWidget_P->graph(0)->setName("рек");
    ui->PlotWidget_P->graph(0)->setPen(pen_advice);

    ui->PlotWidget_P->addGraph();
    ui->PlotWidget_P->graph(1)->setData(t_symm, P_real);
    ui->PlotWidget_P->graph(1)->setName("факт");
    ui->PlotWidget_P->graph(1)->setPen(QPen(Qt::green));

    ui->PlotWidget_P->xAxis->setRange(0, T);
    ui->PlotWidget_P->yAxis->setRange(- 10, P_max + 10);
    ui->PlotWidget_P->xAxis->setLabel("t");
    ui->PlotWidget_P->yAxis->setLabel("P");
    ui->PlotWidget_P->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignBottom);
    ui->PlotWidget_P->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_P->replot();
    ui->PlotWidget_P->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_power_V_explicit.pdf");

    ui->PlotWidget_PT->legend->setVisible(true);

    ui->PlotWidget_PT->addGraph();
    ui->PlotWidget_PT->graph(0)->setData(t_symm, PT_advice);
    ui->PlotWidget_PT->graph(0)->setName("рек");
    ui->PlotWidget_PT->graph(0)->setPen(pen_advice);

    ui->PlotWidget_PT->addGraph();
    ui->PlotWidget_PT->graph(1)->setData(t_symm, P_real);
    ui->PlotWidget_PT->graph(1)->setName("факт");
    ui->PlotWidget_PT->graph(1)->setPen(QPen(Qt::green));

    ui->PlotWidget_PT->xAxis->setRange(0, T);
    ui->PlotWidget_PT->yAxis->setRange(- 10, P_max + 10);
    ui->PlotWidget_PT->xAxis->setLabel("t");
    ui->PlotWidget_PT->yAxis->setLabel("P");
    ui->PlotWidget_PT->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignBottom);
    ui->PlotWidget_PT->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_PT->replot();
    ui->PlotWidget_PT->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_power_T_explicit.pdf");

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

    ui->PlotWidget_N->xAxis->setRange(0, T);
    ui->PlotWidget_N->yAxis->setRange(N_min - (N_max - N_min) * 0.05, N_max + (N_max - N_min) * 0.05);
    ui->PlotWidget_N->xAxis->setLabel("t");
    ui->PlotWidget_N->yAxis->setLabel("N");
    ui->PlotWidget_N->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_N->replot();

    ui->PlotWidget_N->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_N_explicit.pdf");

    double U_max_1 = *std::max_element(U1.begin(), U1.end());
    double U_max_2 = *std::max_element(U2.begin(), U2.end());
    double U_max_3 = *std::max_element(U3.begin(), U3.end());
    double U_min_1 = *std::min_element(U1.begin(), U1.end());
    double U_min_2 = *std::min_element(U2.begin(), U2.end());
    double U_min_3 = *std::min_element(U3.begin(), U3.end());
    double U_max = std::max(std::max(U_max_1, U_max_2), U_max_3);
    double U_min = std::min(std::min(U_min_1, U_min_2), U_min_3);

    ui->PlotWidget_U->legend->setVisible(true);

    ui->PlotWidget_U->addGraph();
    ui->PlotWidget_U->graph(0)->setData(t_symm, U1);
    ui->PlotWidget_U->graph(0)->setName("U_1");
    ui->PlotWidget_U->graph(0)->setPen(QPen(Qt::blue));

    ui->PlotWidget_U->addGraph();
    ui->PlotWidget_U->graph(1)->setData(t_symm, U2);
    ui->PlotWidget_U->graph(1)->setName("U_2");
    ui->PlotWidget_U->graph(1)->setPen(QPen(Qt::magenta));

    ui->PlotWidget_U->addGraph();
    ui->PlotWidget_U->graph(2)->setData(t_symm, U3);
    ui->PlotWidget_U->graph(2)->setName("U_3");
    ui->PlotWidget_U->graph(2)->setPen(QPen(Qt::cyan));

    ui->PlotWidget_U->xAxis->setRange(0, T);
    ui->PlotWidget_U->yAxis->setRange(U_min - (U_max - U_min) * 0.05, U_max + (U_max - U_min) * 0.05);
    ui->PlotWidget_U->xAxis->setLabel("t");
    ui->PlotWidget_U->yAxis->setLabel("U");
    ui->PlotWidget_U->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_U->replot();

    ui->PlotWidget_U->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_U_explicit.pdf");

    double nu_max_1 = *std::max_element(nu_1_symm.begin(), nu_1_symm.end());
    double nu_max_2 = *std::max_element(nu_2_symm.begin(), nu_2_symm.end());
    double nu_min_1 = *std::min_element(nu_1_symm.begin(), nu_1_symm.end());
    double nu_min_2 = *std::min_element(nu_2_symm.begin(), nu_2_symm.end());
    double nu_max_12 = std::max(nu_max_1, nu_max_2);
    double nu_min_12 = std::min(nu_min_1, nu_min_2);


    double nu_min_3 = *std::min_element(nu_3_symm.begin(), nu_3_symm.end());
    double nu_max_3 = *std::max_element(nu_3_symm.begin(), nu_3_symm.end());

    ui->PlotWidget_nu_12->legend->setVisible(true);

    ui->PlotWidget_nu_12->addGraph();
    ui->PlotWidget_nu_12->graph(0)->setData(t_symm, nu_1_symm);
    ui->PlotWidget_nu_12->graph(0)->setName("N_1");
    ui->PlotWidget_nu_12->graph(0)->setPen(QPen(Qt::blue));

    ui->PlotWidget_nu_12->addGraph();
    ui->PlotWidget_nu_12->graph(1)->setData(t_symm, nu_2_symm);
    ui->PlotWidget_nu_12->graph(1)->setName("N_1");
    ui->PlotWidget_nu_12->graph(1)->setPen(QPen(Qt::magenta));

    ui->PlotWidget_nu_12->xAxis->setRange(0, T);
    ui->PlotWidget_nu_12->yAxis->setRange(nu_min_12 - (nu_max_12 - nu_min_12) * 0.05, nu_max_12 + (nu_max_12 - nu_min_12) * 0.05);
    ui->PlotWidget_nu_12->xAxis->setLabel("t");
    ui->PlotWidget_nu_12->yAxis->setLabel("nu");
    ui->PlotWidget_nu_12->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_nu_12->replot();

    ui->PlotWidget_nu_12->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_nu_12_explicit.pdf");

    ui->PlotWidget_nu_3->legend->setVisible(true);

    ui->PlotWidget_nu_3->addGraph();
    ui->PlotWidget_nu_3->graph(0)->setData(t_symm, nu_3_symm);
    ui->PlotWidget_nu_3->graph(0)->setName("N_1");
    ui->PlotWidget_nu_3->graph(0)->setPen(QPen(Qt::cyan));

    ui->PlotWidget_nu_3->xAxis->setRange(0, T);
    ui->PlotWidget_nu_3->yAxis->setRange(nu_min_3 - (nu_max_3 - nu_min_3) * 0.05, nu_max_3 + (nu_max_3 - nu_min_3) * 0.05);
    ui->PlotWidget_nu_3->xAxis->setLabel("t");
    ui->PlotWidget_nu_3->yAxis->setLabel("nu");
    ui->PlotWidget_nu_3->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_nu_3->replot();

    ui->PlotWidget_nu_3->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_nu_3_explicit.pdf");

    double theta_min = *std::min_element(theta_symm.begin(), theta_symm.end());
    double theta_max = *std::max_element(theta_symm.begin(), theta_symm.end());

    ui->PlotWidget_theta->addGraph();
    ui->PlotWidget_theta->graph(0)->setData(t_symm, theta_symm);
    ui->PlotWidget_theta->graph(0)->setName("theta");
    ui->PlotWidget_theta->graph(0)->setPen(QPen(Qt::magenta));

    ui->PlotWidget_theta->xAxis->setRange(0, T);
    ui->PlotWidget_theta->yAxis->setRange(theta_min - (theta_max - theta_min) * 0.05, theta_max + (theta_max - theta_min) * 0.05);
    ui->PlotWidget_theta->xAxis->setLabel("t");
    ui->PlotWidget_theta->yAxis->setLabel("theta");
    ui->PlotWidget_theta->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_theta->replot();

    ui->PlotWidget_theta->savePdf("../custom-omni-wheel-vehicle-control/PICS/reactions_theta_explicit.pdf");
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
