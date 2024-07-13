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
    if (plotted)
    {
        t_noeps.clear();
        nu_1_noeps.clear();
        nu_2_noeps.clear();
        nu_3_noeps.clear();
        x_noeps.clear();
        y_noeps.clear();
        theta_noeps.clear();
        chi_2_noeps.clear();

        t.clear();
        nu_1.clear();
        nu_2.clear();
        nu_3.clear();
        x.clear();
        y.clear();
        theta.clear();
        chi_2.clear();
    }

    T = 120;

    DOPRI8_symmetrical_plot (0, T,
                             t_noeps, nu_1_noeps, nu_2_noeps, nu_3_noeps,
                             x_noeps, y_noeps, theta_noeps, chi_2_noeps, 0);
    DOPRI8_symmetrical_plot (0, T,
                             t, nu_1, nu_2, nu_3,
                             x, y, theta, chi_2, sqrt(4e-5));

    if (plotted)
    {
        ui->PlotWidget_trajectory->clearPlottables();
        ui->PlotWidget_nu_1->clearPlottables();
        ui->PlotWidget_nu_2->clearPlottables();
    }

    trajectory_noeps = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);
    trajectory_eps = new QCPCurve(ui->PlotWidget_trajectory->xAxis, ui->PlotWidget_trajectory->yAxis);

    QVector<QCPCurveData> data_noeps, data_eps;

    QPen pen_noeps(Qt::gray);
    pen_noeps.setStyle(Qt::DashLine);
    QPen pen_eps(Qt::magenta);
    trajectory_noeps->setPen(pen_noeps);
    trajectory_eps->setPen(pen_eps);

    for (int i = 0; i < x_noeps.length(); i++) {
        data_noeps.append(QCPCurveData(i, x_noeps[i], y_noeps[i]));
    }

    for (int i = 0; i < x.length(); i++) {
        data_eps.append(QCPCurveData(i, x[i], y[i]));
    }

    trajectory_noeps->data()->set(data_noeps, true);
    trajectory_eps->data()->set(data_eps, true);

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

    double x_max_1 = *std::max_element(x.begin(), x.end());
    double x_max_2 = *std::max_element(x_noeps.begin(), x_noeps.end());
    double x_max = std::max(x_max_1, x_max_2);
    double x_min_1 = *std::min_element(x.begin(), x.end());
    double x_min_2 = *std::min_element(x_noeps.begin(), x_noeps.end());
    double x_min = std::min(x_min_1, x_min_2);
    double y_max_1 = *std::max_element(y.begin(), y.end());
    double y_max_2 = *std::max_element(y_noeps.begin(), y_noeps.end());
    double y_max = std::max(y_max_1, y_max_2);
    double y_min_1 = *std::min_element(y.begin(), y.end());
    double y_min_2 = *std::min_element(y_noeps.begin(), y_noeps.end());
    double y_min = std::min(y_min_1, y_min_2);


    ui->PlotWidget_trajectory->legend->setVisible(true);
    ui->PlotWidget_trajectory->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignTop);

    ui->PlotWidget_trajectory->xAxis->setRange(x_min - (x_max - x_min) * 0.05, x_max + (x_max - x_min) * 0.05);
    ui->PlotWidget_trajectory->yAxis->setRange(y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05);
    ui->PlotWidget_trajectory->xAxis->setLabel("x");
    ui->PlotWidget_trajectory->yAxis->setLabel("y");

    trajectory_noeps->setName("e^2 = 0 kg*r^2");
    trajectory_eps->setName("e^2 = 4*10^-5 kg*r");

    ui->PlotWidget_trajectory->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_trajectory->replot();

    ui->PlotWidget_trajectory->savePdf("../custom-omni-wheel-vehicle-control/PICS/massive_trajectory.pdf");

    ui->PlotWidget_nu_1->legend->setVisible(true);

    ui->PlotWidget_nu_1->addGraph();
    ui->PlotWidget_nu_1->graph(0)->setData(t_noeps, nu_1_noeps);
    ui->PlotWidget_nu_1->graph(0)->setName("e^2 = 0 kg*r^2");
    ui->PlotWidget_nu_1->graph(0)->setPen(pen_noeps);

    ui->PlotWidget_nu_1->addGraph();
    ui->PlotWidget_nu_1->graph(1)->setData(t, nu_1);
    ui->PlotWidget_nu_1->graph(1)->setName("e^2 = 4*10^-5 kg*r");
    ui->PlotWidget_nu_1->graph(1)->setPen(pen_eps);

    double nu_1_max_1 = *std::max_element(nu_1.begin(), nu_1.end());
    double nu_1_max_2 = *std::max_element(nu_1_noeps.begin(), nu_1_noeps.end());
    double nu_1_max = std::max(nu_1_max_1, nu_1_max_2);
    double nu_1_min_1 = *std::min_element(nu_1.begin(), nu_1.end());
    double nu_1_min_2 = *std::min_element(nu_1_noeps.begin(), nu_1_noeps.end());
    double nu_1_min = std::min(nu_1_min_1, nu_1_min_2);


    ui->PlotWidget_nu_1->xAxis->setRange(0, T);
    ui->PlotWidget_nu_1->yAxis->setRange(nu_1_min - (nu_1_max - nu_1_min) * 0.05, nu_1_max + (nu_1_max - nu_1_min) * 0.05);
    ui->PlotWidget_nu_1->xAxis->setLabel("t");
    ui->PlotWidget_nu_1->yAxis->setLabel("nu_1");
    ui->PlotWidget_nu_1->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignTop);
    ui->PlotWidget_nu_1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_nu_1->replot();
    ui->PlotWidget_nu_1->savePdf("../custom-omni-wheel-vehicle-control/PICS/massive_nu_1.pdf");

    ui->PlotWidget_nu_2->legend->setVisible(true);

    ui->PlotWidget_nu_2->addGraph();
    ui->PlotWidget_nu_2->graph(0)->setData(t_noeps, nu_2_noeps);
    ui->PlotWidget_nu_2->graph(0)->setName("e^2 = 0 kg*r^2");
    ui->PlotWidget_nu_2->graph(0)->setPen(pen_noeps);

    ui->PlotWidget_nu_2->addGraph();
    ui->PlotWidget_nu_2->graph(1)->setData(t, nu_2);
    ui->PlotWidget_nu_2->graph(1)->setName("e^2 = 4*10^-5 kg*r");
    ui->PlotWidget_nu_2->graph(1)->setPen(pen_eps);

    double nu_2_max_1 = *std::max_element(nu_2.begin(), nu_2.end());
    double nu_2_max_2 = *std::max_element(nu_2_noeps.begin(), nu_2_noeps.end());
    double nu_2_max = std::max(nu_2_max_1, nu_2_max_2);
    double nu_2_min_1 = *std::min_element(nu_2.begin(), nu_2.end());
    double nu_2_min_2 = *std::min_element(nu_2_noeps.begin(), nu_2_noeps.end());
    double nu_2_min = std::min(nu_2_min_1, nu_2_min_2);


    ui->PlotWidget_nu_2->xAxis->setRange(0, T);
    ui->PlotWidget_nu_2->yAxis->setRange(nu_2_min - (nu_2_max - nu_2_min) * 0.05, nu_2_max + (nu_2_max - nu_2_min) * 0.05);
    ui->PlotWidget_nu_2->xAxis->setLabel("t");
    ui->PlotWidget_nu_2->yAxis->setLabel("nu_2");
    ui->PlotWidget_nu_2->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignTop);
    ui->PlotWidget_nu_2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_nu_2->replot();
    ui->PlotWidget_nu_2->savePdf("../custom-omni-wheel-vehicle-control/PICS/massive_nu_2.pdf");

    ui->PlotWidget_nu_3->legend->setVisible(true);

    ui->PlotWidget_nu_3->addGraph();
    ui->PlotWidget_nu_3->graph(0)->setData(t_noeps, nu_3_noeps);
    ui->PlotWidget_nu_3->graph(0)->setName("e^2 = 0 kg*r^2");
    ui->PlotWidget_nu_3->graph(0)->setPen(pen_noeps);

    ui->PlotWidget_nu_3->addGraph();
    ui->PlotWidget_nu_3->graph(1)->setData(t, nu_3);
    ui->PlotWidget_nu_3->graph(1)->setName("e^2 = 4*10^-5 kg*r");
    ui->PlotWidget_nu_3->graph(1)->setPen(pen_eps);

    double nu_3_max_1 = *std::max_element(nu_3.begin(), nu_3.end());
    double nu_3_max_2 = *std::max_element(nu_3_noeps.begin(), nu_3_noeps.end());
    double nu_3_max = std::max(nu_3_max_1, nu_3_max_2);
    double nu_3_min_1 = *std::min_element(nu_3.begin(), nu_3.end());
    double nu_3_min_2 = *std::min_element(nu_3_noeps.begin(), nu_3_noeps.end());
    double nu_3_min = std::min(nu_3_min_1, nu_3_min_2);


    ui->PlotWidget_nu_3->xAxis->setRange(0, T);
    ui->PlotWidget_nu_3->yAxis->setRange(nu_3_min - (nu_3_max - nu_3_min) * 0.05, nu_3_max + (nu_3_max - nu_3_min) * 0.05);
    ui->PlotWidget_nu_3->xAxis->setLabel("t");
    ui->PlotWidget_nu_3->yAxis->setLabel("nu_3");
    ui->PlotWidget_nu_3->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignTop);
    ui->PlotWidget_nu_3->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_nu_3->replot();
    ui->PlotWidget_nu_3->savePdf("../custom-omni-wheel-vehicle-control/PICS/massive_nu_3.pdf");

    ui->PlotWidget_chi_2->legend->setVisible(true);

    ui->PlotWidget_chi_2->addGraph();
    ui->PlotWidget_chi_2->graph(0)->setData(t_noeps, chi_2_noeps);
    ui->PlotWidget_chi_2->graph(0)->setName("e^2 = 0 kg*r^2");
    ui->PlotWidget_chi_2->graph(0)->setPen(pen_noeps);

    ui->PlotWidget_chi_2->addGraph();
    ui->PlotWidget_chi_2->graph(1)->setData(t, chi_2);
    ui->PlotWidget_chi_2->graph(1)->setName("e^2 = 4*10^-5 kg*r");
    ui->PlotWidget_chi_2->graph(1)->setPen(pen_eps);

    double chi_2_max_1 = *std::max_element(chi_2.begin(), chi_2.end());
    double chi_2_max_2 = *std::max_element(chi_2_noeps.begin(), chi_2_noeps.end());
    double chi_2_max = std::max(chi_2_max_1, chi_2_max_2);
    double chi_2_min_1 = *std::min_element(chi_2.begin(), chi_2.end());
    double chi_2_min_2 = *std::min_element(chi_2_noeps.begin(), chi_2_noeps.end());
    double chi_2_min = std::min(chi_2_min_1, chi_2_min_2);


    ui->PlotWidget_chi_2->xAxis->setRange(0, T);
    ui->PlotWidget_chi_2->yAxis->setRange(chi_2_min - (chi_2_max - chi_2_min) * 0.05, chi_2_max + (chi_2_max - chi_2_min) * 0.05);
    ui->PlotWidget_chi_2->xAxis->setLabel("t");
    ui->PlotWidget_chi_2->yAxis->setLabel("chi_2");
    ui->PlotWidget_chi_2->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignTop);
    ui->PlotWidget_chi_2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->PlotWidget_chi_2->replot();
    ui->PlotWidget_chi_2->savePdf("../custom-omni-wheel-vehicle-control/PICS/massive_chi_2.pdf");

    plotted = true;
}
