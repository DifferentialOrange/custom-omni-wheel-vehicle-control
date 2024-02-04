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

const double program_R = 2;

Vector<2> get_program_nu3_coeff(double theta_square_coef, double theta_lin_coef) {
    Vector<2> res;

    res[0] = 2 * M_PI * parameters::symmetrical::L * theta_square_coef;
    res[1] = 2 * M_PI * parameters::symmetrical::L * theta_lin_coef;
//    res[0] = 0;
//    res[1] = 2 * M_PI * parameters::symmetrical::L * thetaN_1;

    qDebug() << "nu3(t) = " << res[0] << "t +" << res[1] << "\n";

    return res;
}

Vector<3> get_expected_nu(double t, double theta_square_coef, double theta_lin_coef) {
    Vector<3> res;

    auto nu_coeff = get_program_nu3_coeff(theta_square_coef, theta_lin_coef);

    res[0] = 0;
    res[2] = nu_coeff[0] * t + nu_coeff[1];
    res[1] = program_R * res[2] / parameters::symmetrical::L;

    return res;
}

Vector<3> get_expected_pos(double t, double theta_square_coef, double theta_lin_coef) {
    auto nu_coeff = get_program_nu3_coeff(theta_square_coef, theta_lin_coef);

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
            qDebug() << "there were " << loop_iter(theta[i - 1]) + 1 << "completed loops\n";
            return i;
        }
    }

    return 0;
}

void MainWindow::on_pushButton_compute_clicked()
{
    bool ok;

    double nu_3_lin_coef = 0.01;
    if (!ok)
        return;

    double nu_3_const_coef = 0.1;
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

        U1.clear();
        U2.clear();
        U3.clear();
    int steps_per_loop = 30;
    double theta_square_coef = nu_3_lin_coef;
    double theta_lin_coef = nu_3_const_coef;
    double first_loop_length = (
                (sqrt(theta_lin_coef * theta_lin_coef + 2 * theta_square_coef) - theta_lin_coef)
            / // -------------------------------------------------------------------------------
                                          theta_square_coef
                );
    qDebug() << "sqrt(theta_lin_coef * theta_lin_coef + 2 * theta_square_coef) " << sqrt(theta_lin_coef * theta_lin_coef + 2 * theta_square_coef) << "\n";
    qDebug() << "theta_lin_coef " << theta_lin_coef << "\n";
    qDebug() << "theta_lin_coef " << theta_lin_coef << "\n";
    qDebug() << "first_loop_length " << first_loop_length << "\n";
    int steps = steps_per_loop * (T / first_loop_length);

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

            auto expected_start_nu = get_expected_nu(real_initial_t, theta_square_coef, theta_lin_coef);
            auto expected_start_pos = get_expected_pos(real_initial_t, theta_square_coef, theta_lin_coef);

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

        auto expected_nu = get_expected_nu(real_final_t, theta_square_coef, theta_lin_coef);
        auto expected_pos = get_expected_pos(real_final_t, theta_square_coef, theta_lin_coef);

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

        if (!model_satisfied(N_1, N_2, N_3)) {
            break;
        }
    }

    DOPRI8_symmetrical_plot (0, T, initial_values, control_1,
                             control_2, t_sw,
                             t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                             x_symm, y_symm, theta_symm,
                             P_real, P_advice, N_1, N_2, N_3, U1, U2, U3);

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_P->clearPlottables();
        ui->PlotWidget_N->clearPlottables();
        ui->PlotWidget_U->clearPlottables();
        ui->PlotWidget_nu->clearPlottables();
    }

    trajectory_minus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_minus_symm, data_plus_symm;

    QPen pen_minus_symm(Qt::gray);
    QPen pen_plus_symm(Qt::gray);
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

    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/trajectory_explicit.pdf");

    double P_max_1 = *std::max_element(P_real.begin(), P_real.end());
    double P_max_2 = *std::max_element(P_advice.begin(), P_advice.end());
    double P_max = std::max(P_max_1, P_max_2);

    ui->PlotWidget_P->legend->setVisible(true);

    QPen pen_advice(Qt::red);
    pen_advice.setStyle(Qt::DashLine);
    ui->PlotWidget_P->addGraph();
    ui->PlotWidget_P->graph(0)->setData(t_symm, P_advice);
    ui->PlotWidget_P->graph(0)->setName("advice");
    ui->PlotWidget_P->graph(0)->setPen(pen_advice);

    ui->PlotWidget_P->addGraph();
    ui->PlotWidget_P->graph(1)->setData(t_symm, P_real);
    ui->PlotWidget_P->graph(1)->setName("real");
    ui->PlotWidget_P->graph(1)->setPen(QPen(Qt::green));

    ui->PlotWidget_P->xAxis->setRange(0, T);
    ui->PlotWidget_P->yAxis->setRange(- 0.05, P_max + 0.05);
    ui->PlotWidget_P->xAxis->setLabel("t");
    ui->PlotWidget_P->yAxis->setLabel("P");
    ui->PlotWidget_P->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_P->replot();
    ui->PlotWidget_P->savePdf("../custom-omni-wheel-vehicle-control/PICS/power_explicit.pdf");

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
    ui->PlotWidget_N->yAxis->setRange(N_min - 0.05, N_max + 0.05);
    ui->PlotWidget_N->xAxis->setLabel("t");
    ui->PlotWidget_N->yAxis->setLabel("N");
    ui->PlotWidget_N->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_N->replot();

    ui->PlotWidget_N->savePdf("../custom-omni-wheel-vehicle-control/PICS/N_explicit.pdf");

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
    ui->PlotWidget_U->yAxis->setRange(U_min - 0.05, U_max + 0.05);
    ui->PlotWidget_U->xAxis->setLabel("t");
    ui->PlotWidget_U->yAxis->setLabel("U");
    ui->PlotWidget_U->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_U->replot();

    ui->PlotWidget_U->savePdf("../custom-omni-wheel-vehicle-control/PICS/U_explicit.pdf");

    double nu_max_1 = *std::max_element(nu_1_symm.begin(), nu_1_symm.end());
    double nu_max_2 = *std::max_element(nu_2_symm.begin(), nu_2_symm.end());
    double nu_max_3 = *std::max_element(nu_3_symm.begin(), nu_3_symm.end());
    double nu_min_1 = *std::min_element(nu_1_symm.begin(), nu_1_symm.end());
    double nu_min_2 = *std::min_element(nu_2_symm.begin(), nu_2_symm.end());
    double nu_min_3 = *std::min_element(nu_3_symm.begin(), nu_3_symm.end());
    double nu_max = std::max(std::max(nu_max_1, nu_max_2), nu_max_3);
    double nu_min = std::min(std::min(nu_min_1, nu_min_2), nu_min_3);

    ui->PlotWidget_nu->legend->setVisible(true);

    ui->PlotWidget_nu->addGraph();
    ui->PlotWidget_nu->graph(0)->setData(t_symm, nu_1_symm);
    ui->PlotWidget_nu->graph(0)->setName("nu_1");
    ui->PlotWidget_nu->graph(0)->setPen(QPen(Qt::blue));

    ui->PlotWidget_nu->addGraph();
    ui->PlotWidget_nu->graph(1)->setData(t_symm, nu_2_symm);
    ui->PlotWidget_nu->graph(1)->setName("nu_2");
    ui->PlotWidget_nu->graph(1)->setPen(QPen(Qt::magenta));

    ui->PlotWidget_nu->addGraph();
    ui->PlotWidget_nu->graph(2)->setData(t_symm, nu_3_symm);
    ui->PlotWidget_nu->graph(2)->setName("nu_3");
    ui->PlotWidget_nu->graph(2)->setPen(QPen(Qt::cyan));

    ui->PlotWidget_nu->xAxis->setRange(0, T);
    ui->PlotWidget_nu->yAxis->setRange(nu_min - 0.05, nu_max + 0.05);
    ui->PlotWidget_nu->xAxis->setLabel("t");
    ui->PlotWidget_nu->yAxis->setLabel("nu");
    ui->PlotWidget_nu->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_nu->replot();

    ui->PlotWidget_nu->savePdf("../custom-omni-wheel-vehicle-control/PICS/nu_explicit.pdf");

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
