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
    ui->lineEdit_x_T->setText("5");
    ui->lineEdit_y_T->setText("5");
    ui->lineEdit_theta_T->setText("0");
    ui->lineEdit_t_sw->setText("5");
    ui->lineEdit_T->setText("10");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void fillTrajectory(double t_sw, QVector<double> &t, QVector<double> &t_symm,
                    QVector<double> &x, QVector<double> &x_symm,
                    QVector<double> &y, QVector<double> &y_symm,
                    QCustomPlot* window) {

    QCPCurve* trajectory_minus_symm = new QCPCurve(window->xAxis, window->yAxis);
    QCPCurve* trajectory_plus_symm = new QCPCurve(window->xAxis, window->yAxis);
    QCPCurve* trajectory_minus = new QCPCurve(window->xAxis, window->yAxis);
    QCPCurve* trajectory_plus = new QCPCurve(window->xAxis, window->yAxis);

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

    trajectory_minus->setName("e^2");
    trajectory_plus->setName("e^2");
}

void fillGraph(double t_sw, QVector<double> &t, QVector<double> &t_symm,
                    QVector<double> &v, QVector<double> &v_symm,
                    QCustomPlot* window) {
    QCPCurve* v_minus_symm = new QCPCurve(window->xAxis, window->yAxis);
    QCPCurve* v_plus_symm = new QCPCurve(window->xAxis, window->yAxis);
    QCPCurve* v_minus = new QCPCurve(window->xAxis, window->yAxis);
    QCPCurve* v_plus = new QCPCurve(window->xAxis, window->yAxis);

    QVector<QCPCurveData> data_minus, data_plus, data_minus_symm, data_plus_symm;

    QPen pen_minus_symm(Qt::DashLine);
    pen_minus_symm.setColor(Qt::gray);
    QPen pen_plus_symm(Qt::DashLine);
    pen_plus_symm.setColor(Qt::yellow);
    v_minus_symm->setPen(pen_minus_symm);
    v_plus_symm->setPen(pen_plus_symm);

    QPen pen_minus(Qt::blue);
    QPen pen_plus(Qt::magenta);
    v_minus->setPen(pen_minus);
    v_plus->setPen(pen_plus);

    int i = 0;
    for (i = 0; t_symm[i] < t_sw; i++)
        data_minus_symm.append(QCPCurveData(i, t_symm[i], v_symm[i]));

    // To plot voltages without gap.
    data_plus_symm.append(QCPCurveData(i, t_symm[i - 1], v_symm[i - 1]));

    for (; i < t_symm.length(); i++)
        data_plus_symm.append(QCPCurveData(i, t_symm[i], v_symm[i]));

    v_minus_symm->data()->set(data_minus_symm, true);
    v_plus_symm->data()->set(data_plus_symm, true);

    for (i = 0; t[i] < t_sw; i++)
        data_minus.append(QCPCurveData(i, t[i], v[i]));

    // To plot voltages without gap.
    data_plus.append(QCPCurveData(i, t[i - 1], v[i - 1]));

    for (; i < t.length(); i++)
        data_plus.append(QCPCurveData(i, t[i], v[i]));

    v_minus->data()->set(data_minus, true);
    v_plus->data()->set(data_plus, true);

    v_minus->setName("e^2");
    v_plus->setName("e^2");
}

void setRange(QCPAxis *ax, QVector<double> &v_1, QVector<double> &v_2, double margin) {
    double v_max_1 = *std::max_element(v_1.begin(), v_1.end());
    double v_min_1 = *std::min_element(v_1.begin(), v_1.end());

    double v_max_2 = *std::max_element(v_2.begin(), v_2.end());
    double v_min_2 = *std::min_element(v_2.begin(), v_2.end());

    double v_max = std::max(v_max_1, v_max_2);
    double v_min = std::min(v_min_1, v_min_2);

    ax->setRange(v_min - (v_max - v_min) * margin, v_max + (v_max - v_min) * margin);
}

void doGraph(double t_sw, QVector<double> &t, QVector<double> &t_symm,
             QVector<double> &v, QVector<double> &v_symm,
             QCustomPlot* window, QString filename,
             double T, Vector<6> initial_values, Vector<6> final_values) {
    fillGraph(t_sw, t, t_symm, v, v_symm, window);
    setRange(window->xAxis, t_symm, t, 0);
    setRange(window->yAxis, v_symm, v, 0.05);

    window->xAxis->setLabel("t");
    window->yAxis->setLabel(filename);

    window->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    window->replot();

    window->savePdf("../custom-omni-wheel-vehicle-control/PICS/" + filename + "_t_sw_"
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
        u_1_symm.clear();
        u_2_symm.clear();
        u_3_symm.clear();

        t.clear();
        nu_1.clear();
        nu_2.clear();
        nu_3.clear();
        x.clear();
        y.clear();
        theta.clear();
        u_1.clear();
        u_2.clear();
        u_3.clear();
    }

    Vector<6> control = predict_control(t_sw, T, initial_values[0], final_values[0],
            initial_values[1], final_values[1], initial_values[2], final_values[2],
            final_values[3], final_values[4], final_values[5]);

    Vector<6> u_symm = control;

    DOPRI8_symmetrical_plot (0, T, initial_values, {control[0], control[1], control[2]},
                             {control[3], control[4], control[5]}, t_sw,
                             t_symm, nu_1_symm, nu_2_symm, nu_3_symm,
                             x_symm, y_symm, theta_symm,
                             u_1_symm, u_2_symm, u_3_symm);

    control = custom_control_find(control, t_sw, T, initial_values, final_values);

    DOPRI8_final_plot (0, T, initial_values, {control[0], control[1], control[2]},
                        {control[3], control[4], control[5]},
                        t_sw, t, nu_1, nu_2, nu_3, x, y, theta, u_1, u_2, u_3);

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_theta->clearPlottables();
        ui->PlotWidget_nu_1->clearPlottables();
        ui->PlotWidget_nu_2->clearPlottables();
        ui->PlotWidget_nu_3->clearPlottables();
        ui->PlotWidget_u_1->clearPlottables();
        ui->PlotWidget_u_2->clearPlottables();
        ui->PlotWidget_u_3->clearPlottables();
    }

    fillTrajectory(t_sw, t, t_symm, x, x_symm, y, y_symm, ui->PlotWidget_trajectory);
    setRange(ui->PlotWidget_trajectory->xAxis, x_symm, x, 0.05);
    setRange(ui->PlotWidget_trajectory->yAxis, y_symm, y, 0.05);

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

    fillGraph(t_sw, t, t_symm, theta, theta_symm, ui->PlotWidget_theta);
    setRange(ui->PlotWidget_theta->xAxis, t_symm, t, 0);
    setRange(ui->PlotWidget_theta->yAxis, theta_symm, theta, 0.05);

    ui->PlotWidget_theta->xAxis->setLabel("t");
    ui->PlotWidget_theta->yAxis->setLabel("th");

    ui->PlotWidget_theta->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_theta->replot();

    ui->PlotWidget_theta->savePdf("../custom-omni-wheel-vehicle-control/PICS/theta_t_sw_"
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

    doGraph(t_sw, t, t_symm, theta, theta_symm, ui->PlotWidget_theta, "theta", T, initial_values, final_values);
    doGraph(t_sw, t, t_symm, nu_1, nu_1_symm, ui->PlotWidget_nu_1, "nu_1", T, initial_values, final_values);
    doGraph(t_sw, t, t_symm, nu_2, nu_2_symm, ui->PlotWidget_nu_2, "nu_2", T, initial_values, final_values);
    doGraph(t_sw, t, t_symm, nu_3, nu_3_symm, ui->PlotWidget_nu_3, "nu_3", T, initial_values, final_values);
    doGraph(t_sw, t, t_symm, u_1, u_1_symm, ui->PlotWidget_u_1, "u_1", T, initial_values, final_values);
    doGraph(t_sw, t, t_symm, u_2, u_2_symm, ui->PlotWidget_u_2, "u_2", T, initial_values, final_values);
    doGraph(t_sw, t, t_symm, u_3, u_3_symm, ui->PlotWidget_u_3, "u_3", T, initial_values, final_values);

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
