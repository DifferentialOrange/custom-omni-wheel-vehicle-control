#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include.h"
#include <fstream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
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

static double energy(Vector<6> control, double t_sw, double T)
{
    return t_sw * (control[0] * control[0] + control[1] * control[1] + control[2] * control[2]) + \
            (T - t_sw) * (control[3] * control[3] + control[4] * control[4] + control[5] * control[5]);
}

void MainWindow::on_pushButton_compute_clicked()
{
    std::ofstream energy_0("energy_0.dat");
    std::ofstream energy_turn("energy_turn.dat");

    initial_values[0] = 0;
    initial_values[1] = 0;
    initial_values[2] = 0;
    initial_values[3] = 0;
    initial_values[4] = 0;
    initial_values[5] = 0;
    final_values[0] = 0;
    final_values[1] = 0;
    final_values[2] = 0;

    double T = 1;

    for (final_values[3] = -10; final_values[3] < 10.1; final_values[3] += 0.1)
         for (final_values[4] = -10; final_values[4] < 10.1; final_values[4] += 0.1)
         {
             double en_0_min = 1e30;
             double en_turn_min = 1e30;

             final_values[5] = 0;

             for (double t_sw = T / 20; t_sw < (T * 99 / 100); t_sw += T / 20)
             {
                 auto control = predict_control(t_sw, T, initial_values[0], final_values[0], initial_values[1], final_values[1],
                         initial_values[2], final_values[2], final_values[3], final_values[4], final_values[5]);

                 en_0_min = std::min(en_0_min, energy(control, t_sw, T));
             }

             energy_0 << final_values[3] << " " << final_values[4] << " " << en_0_min << std::endl;

             final_values[5] = atan2(final_values[4], final_values[3]);

             for (double t_sw = T / 20; t_sw < (T * 99 / 100); t_sw += T / 20)
             {
                 auto control = predict_control(t_sw, T, initial_values[0], final_values[0], initial_values[1], final_values[1],
                         initial_values[2], final_values[2], final_values[3], final_values[4], final_values[5]);

                 en_turn_min = std::min(en_turn_min, energy(control, t_sw, T));
             }

             energy_turn << final_values[3] << " " << final_values[4] << " " << en_turn_min << std::endl;

             qDebug() << "Finished (" << final_values[3] << ", " << final_values[4] << ")\n";
         }
}
