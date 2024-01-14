#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    plotted(false)
{
    ui->setupUi(this);

    ui->lineEdit_N_1->setText("0.1");
    ui->lineEdit_N_2->setText("0.21");
    ui->lineEdit_T->setText("20");

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

const double program_R = 2;

Vector<2> get_program_nu3_coeff(double thetaN_1, double thetaN_2) {
    Vector<2> res;

    res[0] = 2 * M_PI * parameters::symmetrical::L * (thetaN_2 - 2 * thetaN_1);
    res[1] = 2 * M_PI * parameters::symmetrical::L * (2 * thetaN_1 - thetaN_2 / 2);
//    res[0] = 0;
//    res[1] = 2 * M_PI * parameters::symmetrical::L * thetaN_1;

    return res;
}

Vector<3> get_expected_nu(double t, double thetaN_1, double thetaN_2) {
    Vector<3> res;

    auto nu_coeff = get_program_nu3_coeff(thetaN_1, thetaN_2);

    res[0] = 0;
    res[2] = nu_coeff[0] * t + nu_coeff[1];
    res[1] = program_R * res[2] / parameters::symmetrical::L;

    return res;
}

Vector<3> get_expected_pos(double t, double thetaN_1, double thetaN_2) {
    auto nu_coeff = get_program_nu3_coeff(thetaN_1, thetaN_2);

    double theta = (nu_coeff[0] * t * t / 2 + nu_coeff[1] * t) / parameters::symmetrical::L;
    double x = program_R * (cos(theta) - 1);
    double y = program_R * sin(theta);

    Vector<3> res;

    res[0] = x;
    res[1] = y;
    res[2] = theta;

    return res;
}

bool N_satisfies_the_model(QVector<double> N) {
    return *std::min_element(N.begin(), N.end()) > 0;
}

bool model_satisfied(QVector<double> N_1, QVector<double> N_2, QVector<double> N_3) {
    return N_satisfies_the_model(N_1) && N_satisfies_the_model(N_2) && N_satisfies_the_model(N_3);
}

int div_ceil(double d1, double d2) {
    return ceil(d1 / d2);
}

double loop_iter(double theta) {
    return div_ceil(theta, 2 * M_PI);
}

int find_last_loop_start(QVector<double> theta) {
    for (int i = theta.length() - 1; i > 0; i--) {
        if (loop_iter(theta[i]) != loop_iter(theta[i - 1])) {
            return i;
        }
    }

    return 0;
}

void MainWindow::on_pushButton_compute_clicked()
{
    bool ok;

    double thetaN_1 = ui->lineEdit_N_1->text().toDouble(&ok);
    if (!ok)
        return;

    double thetaN_2 = ui->lineEdit_N_2->text().toDouble(&ok);
    if (!ok)
        return;

    double T = ui->lineEdit_T->text().toDouble(&ok);
    if (!ok)
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

    int steps_per_loop = 30;
    double first_loop_length = 1 / thetaN_1;
    int steps = steps_per_loop * (T / first_loop_length);

    bool break_next = false;

    for (int i = 0; i < steps; i++) {
        qDebug() << "iteration " << i << "\n";
        double t_step = T / steps;

        double real_initial_t;
        double real_t_sw;
        double real_final_t;
        Vector<6> real_initial_values;
        Vector<6> real_final_values;

        if (i == 0) {
            real_initial_t = 0.0;

            auto expected_start_nu = get_expected_nu(real_initial_t, thetaN_1, thetaN_2);
            auto expected_start_pos = get_expected_pos(real_initial_t, thetaN_1, thetaN_2);

            real_initial_values[0] = expected_start_nu[0];
            real_initial_values[1] = expected_start_nu[1];
            real_initial_values[2] = expected_start_nu[2];
            real_initial_values[3] = expected_start_pos[0];
            real_initial_values[4] = expected_start_pos[1];
            real_initial_values[5] = expected_start_pos[2];
        } else {
            real_initial_t = t_symm.last();

            real_initial_values[0] = nu_1_symm.last();
            real_initial_values[1] = nu_2_symm.last();
            real_initial_values[2] = nu_3_symm.last();
            real_initial_values[3] = x_symm.last();
            real_initial_values[4] = y_symm.last();
            real_initial_values[5] = theta_symm.last();
        }

        real_t_sw = real_initial_t + t_step / 2;
        real_final_t = real_initial_t + t_step;

        auto expected_nu = get_expected_nu(real_final_t, thetaN_1, thetaN_2);
        auto expected_pos = get_expected_pos(real_final_t, thetaN_1, thetaN_2);

        real_final_values[0] = expected_nu[0];
        real_final_values[1] = expected_nu[1];
        real_final_values[2] = expected_nu[2];
        real_final_values[3] = expected_pos[0];
        real_final_values[4] = expected_pos[1];
        real_final_values[5] = expected_pos[2];

        qDebug() << "from " << "x: " << real_initial_values[3] << "y: " << real_initial_values[4] << "theta: " << real_initial_values[5] << "\n";
        qDebug() << "to " << "x: " << real_final_values[3] << "y: " << real_final_values[4] << "theta: " << real_final_values[5] << "\n";

        double relative_initial_t = 0;
        double relative_t_sw = real_t_sw - real_initial_t;
        double relative_final_t = real_final_t - real_initial_t;

        Vector<6> relative_initial_values;
        relative_initial_values[0] = real_initial_values[0],
        relative_initial_values[1] = real_initial_values[1];
        relative_initial_values[2] = real_initial_values[2];
        relative_initial_values[3] = 0;
        relative_initial_values[4] = 0;
//        relative_initial_values[5] = 0;
        relative_initial_values[5] = real_initial_values[5];

        Vector<6> relative_final_values;
        relative_final_values[0] = real_final_values[0];
        relative_final_values[1] = real_final_values[1];
        relative_final_values[2] = real_final_values[2];
        relative_final_values[3] = real_final_values[3] - real_initial_values[3];
        relative_final_values[4] = real_final_values[4] - real_initial_values[4];
//        relative_final_values[5] = real_final_values[5] - real_initial_values[5];
        relative_final_values[5] = real_final_values[5];

        Vector<6> control = predict_control(
                    relative_t_sw,
                    relative_final_t,
                    relative_initial_values[0],
                    relative_final_values[0],
                    relative_initial_values[1],
                    relative_final_values[1],
                    relative_initial_values[2],
                    relative_final_values[2],
                    relative_final_values[3],
                    relative_final_values[4],
                    relative_initial_values[5],
                    relative_final_values[5]);

        Vector<3> control_1 = {control[0], control[1], control[2]};
        Vector<3> control_2 = {control[3], control[4], control[5]};

        DOPRI8_symmetrical_plot (real_initial_t,
                                 real_final_t,
                                 real_initial_values,
                                 control_1,
                                 control_2,
                                 real_t_sw,
                                 t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                                 x_symm, y_symm, theta_symm,
                                 P_real, P_advice, N_1, N_2, N_3);

        if (break_next) {
            break;
        }

        if (!model_satisfied(N_1, N_2, N_3)) {
            qDebug() << "break next" << "\n";
            break_next = true;
        }
    }

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_P->clearPlottables();
        ui->PlotWidget_N->clearPlottables();
    }

    // Use different color for the last loop;
    trajectory_minus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_minus_symm, data_plus_symm;

    QPen pen_minus_symm(Qt::blue);
    QPen pen_plus_symm(Qt::magenta);
    pen_plus_symm.setWidth(pen_plus_symm.width() + 1);
    trajectory_minus_symm->setPen(pen_minus_symm);
    trajectory_plus_symm->setPen(pen_plus_symm);

    int i = 0;
    int i_boundary = find_last_loop_start(theta_symm);

    qDebug() << "last loop start " << theta_symm[i_boundary] << "last loop end " << theta_symm.last() << "\n";

    for (i = 0; i < i_boundary; i++)
        data_minus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

    for (; i < x_symm.length(); i++)
        data_plus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

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
    qDebug() << "PlotWidget_N" << "\n";

    plotted = true;
}
