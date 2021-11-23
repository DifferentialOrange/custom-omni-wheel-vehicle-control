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

    parameters::mu_n = ui->lineEdit_mu_n->text().toDouble(&ok);
    if (!ok /*|| parameters::mu_n <= 0*/)
        return;

    parameters::mu_tau = ui->lineEdit_mu_tau->text().toDouble(&ok);
    if (!ok || parameters::mu_tau <= 0)
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

        t.clear();
        nu_1.clear();
        nu_2.clear();
        nu_3.clear();
        x.clear();
        y.clear();
        theta.clear();

        v_sign_n_1.clear();
        v_sign_n_2.clear();
        v_sign_n_3.clear();
        v_sign_tau_1.clear();
        v_sign_tau_2.clear();
        v_sign_tau_3.clear();

        N_1.clear();
        N_2.clear();
        N_3.clear();
    }

    Vector<6> control = predict_control(t_sw, T, initial_values[0], final_values[0],
            initial_values[1], final_values[1], initial_values[2], final_values[2],
            final_values[3], final_values[4], final_values[5]);

    qDebug() << "Control found: " << control[0] << ' ' << control[1] << ' '
             << control[2] << ' ' << control[3] << ' '
             << control[4] << ' ' << control[5] << '\n';

    DOPRI8_symmetrical_plot (0, T, initial_values, {control[0], control[1], control[2]},
                             {control[3], control[4], control[5]}, t_sw,
                             t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                             x_symm, y_symm, theta_symm);

    qDebug() << "Symm plot ok" << '\n';

    DOPRI8_friction_plot (0, T, initial_values, {control[0], control[1], control[2]},
                            {control[3], control[4], control[5]}, t_sw,
                            t, nu_1, nu_2, nu_3,
                            x, y, theta, d_chi_1, d_chi_2, d_chi_3, d_theta,
                            v_sign_tau_1, v_sign_tau_2, v_sign_tau_3,
                            v_sign_n_1, v_sign_n_2, v_sign_n_3,
                            N_1, N_2, N_3);
    qDebug() << "Friction plot ok" << '\n';

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_speed_coef_tau->clearPlottables();
        ui->PlotWidget_speed_coef_n->clearPlottables();
        ui->PlotWidget_N->clearPlottables();
    }

//    TO DO: show error in text panel

    trajectory_minus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus_symm = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_minus = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_plus = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_minus, data_plus, data_minus_symm, data_plus_symm;

    QPen pen_minus_symm(Qt::DashLine);
    pen_minus_symm.setColor(Qt::gray);
    QPen pen_plus_symm(Qt::DashLine);
    pen_plus_symm.setColor(Qt::yellow);
    trajectory_minus_symm->setPen(pen_minus_symm);
    trajectory_plus_symm->setPen(pen_plus_symm);

    QPen pen_minus(Qt::blue);
    QPen pen_plus(Qt::magenta);
    trajectory_minus->setPen(pen_minus);
    trajectory_plus->setPen(pen_plus);

    int i = 0;
    for (i = 0; t_symm[i] < t_sw; i++)
        data_minus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

    for (; i < x_symm.length(); i++)
        data_plus_symm.append(QCPCurveData(i, x_symm[i], y_symm[i]));

    trajectory_minus_symm->data()->set(data_minus_symm, true);
    trajectory_plus_symm->data()->set(data_plus_symm, true);

    for (i = 0; t[i] < t_sw; i++)
        data_minus.append(QCPCurveData(i, x[i], y[i]));

    for (; i < x.length(); i++)
        data_plus.append(QCPCurveData(i, x[i], y[i]));

    trajectory_minus->data()->set(data_minus, true);
    trajectory_plus->data()->set(data_plus, true);

    double x_max_1 = *std::max_element(x_symm.begin(), x_symm.end());
    double x_min_1 = *std::min_element(x_symm.begin(), x_symm.end());
    double y_max_1 = *std::max_element(y_symm.begin(), y_symm.end());
    double y_min_1 = *std::min_element(y_symm.begin(), y_symm.end());

    double x_max_2 = *std::max_element(x.begin(), x.end());
    double x_min_2 = *std::min_element(x.begin(), x.end());
    double y_max_2 = *std::max_element(y.begin(), y.end());
    double y_min_2 = *std::min_element(y.begin(), y.end());

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

    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/trajectory_friction_t_sw_"
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
                                + "_mu_n_" + QString::number(parameters::mu_n, 'g', 6)
                                + "_mu_tau_" + QString::number(parameters::mu_tau, 'g', 6) + ".pdf");


    ui->PlotWidget_speed_coef_tau->legend->setVisible(true);

    ui->PlotWidget_speed_coef_tau->addGraph();
    ui->PlotWidget_speed_coef_tau->graph(0)->setData(t, v_sign_tau_1);
    ui->PlotWidget_speed_coef_tau->graph(0)->setName("v_1");
    ui->PlotWidget_speed_coef_tau->graph(0)->setPen(QPen(Qt::red));

    ui->PlotWidget_speed_coef_tau->addGraph();
    ui->PlotWidget_speed_coef_tau->graph(1)->setData(t, v_sign_tau_2);
    ui->PlotWidget_speed_coef_tau->graph(1)->setName("v_2");
    ui->PlotWidget_speed_coef_tau->graph(1)->setPen(QPen(Qt::green));

    ui->PlotWidget_speed_coef_tau->addGraph();
    ui->PlotWidget_speed_coef_tau->graph(2)->setData(t, v_sign_tau_3);
    ui->PlotWidget_speed_coef_tau->graph(2)->setName("v_3");
    ui->PlotWidget_speed_coef_tau->graph(2)->setPen(QPen(Qt::blue));

    ui->PlotWidget_speed_coef_tau->xAxis->setRange(0, T);
    ui->PlotWidget_speed_coef_tau->yAxis->setRange(-1.05, 1.05);
    ui->PlotWidget_speed_coef_tau->xAxis->setLabel("t");
    ui->PlotWidget_speed_coef_tau->yAxis->setLabel("v_tau_sign");
    ui->PlotWidget_speed_coef_tau->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_speed_coef_tau->replot();


    ui->PlotWidget_speed_coef_n->legend->setVisible(true);

    ui->PlotWidget_speed_coef_n->addGraph();
    ui->PlotWidget_speed_coef_n->graph(0)->setData(t, v_sign_n_1);
    ui->PlotWidget_speed_coef_n->graph(0)->setName("v_1");
    ui->PlotWidget_speed_coef_n->graph(0)->setPen(QPen(Qt::red));

    ui->PlotWidget_speed_coef_n->addGraph();
    ui->PlotWidget_speed_coef_n->graph(1)->setData(t, v_sign_n_2);
    ui->PlotWidget_speed_coef_n->graph(1)->setName("v_2");
    ui->PlotWidget_speed_coef_n->graph(1)->setPen(QPen(Qt::green));

    ui->PlotWidget_speed_coef_n->addGraph();
    ui->PlotWidget_speed_coef_n->graph(2)->setData(t, v_sign_n_3);
    ui->PlotWidget_speed_coef_n->graph(2)->setName("v_3");
    ui->PlotWidget_speed_coef_n->graph(2)->setPen(QPen(Qt::blue));

    ui->PlotWidget_speed_coef_n->xAxis->setRange(0, T);
    ui->PlotWidget_speed_coef_n->yAxis->setRange(-1.05, 1.05);
    ui->PlotWidget_speed_coef_n->xAxis->setLabel("t");
    ui->PlotWidget_speed_coef_n->yAxis->setLabel("v_n_sign");
    ui->PlotWidget_speed_coef_n->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_speed_coef_n->replot();

    double N_1_max = *std::max_element(N_1.begin(), N_1.end());
    double N_2_max = *std::max_element(N_2.begin(), N_2.end());
    double N_3_max = *std::max_element(N_3.begin(), N_3.end());
    double N_max = std::max(N_1_max, std::max(N_2_max, N_3_max));

    double N_1_min = *std::min_element(N_1.begin(), N_1.end());
    double N_2_min = *std::min_element(N_2.begin(), N_2.end());
    double N_3_min = *std::min_element(N_3.begin(), N_3.end());
    double N_min = std::min(N_1_min, std::min(N_2_min, N_3_min));

    ui->PlotWidget_N->legend->setVisible(true);

    ui->PlotWidget_N->addGraph();
    ui->PlotWidget_N->graph(0)->setData(t, N_1);
    ui->PlotWidget_N->graph(0)->setName("N_1");
    ui->PlotWidget_N->graph(0)->setPen(QPen(Qt::red));

    ui->PlotWidget_N->addGraph();
    ui->PlotWidget_N->graph(1)->setData(t, N_2);
    ui->PlotWidget_N->graph(1)->setName("N_2");
    ui->PlotWidget_N->graph(1)->setPen(QPen(Qt::green));

    ui->PlotWidget_N->addGraph();
    ui->PlotWidget_N->graph(2)->setData(t, N_3);
    ui->PlotWidget_N->graph(2)->setName("N_3");
    ui->PlotWidget_N->graph(2)->setPen(QPen(Qt::blue));

    ui->PlotWidget_N->xAxis->setRange(0, T);
    ui->PlotWidget_N->yAxis->setRange(N_min - 0.05, N_max + 0.05);
    ui->PlotWidget_N->xAxis->setLabel("t");
    ui->PlotWidget_N->yAxis->setLabel("N");
    ui->PlotWidget_N->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_N->replot();

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
