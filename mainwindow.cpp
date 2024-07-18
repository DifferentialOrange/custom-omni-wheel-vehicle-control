#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "include.h"
#include <fstream>

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

    std::ofstream file(("max_U_theta_T_"+ QString::number(final_values[5], 'g', 4) + ".txt").toStdString());
    file << "Writing this to a file.\n";

    for (double T = 1.0; T <= 50.0; T += 1.0) {
        for (double t_sw = 0.5; t_sw < T; t_sw += 1.0) {
            Vector<6> control = predict_control(t_sw, T, initial_values[0], final_values[0],
                    initial_values[1], final_values[1], initial_values[2], final_values[2],
                    final_values[3], final_values[4], final_values[5]);

            double u = std::max({fabs(control[0]), fabs(control[1]), fabs(control[2]), fabs(control[3]), fabs(control[4]), fabs(control[5])});
            file << t_sw << " " << T << " " << u << "\n";
        }
    }
    Vector<6> control = predict_control(t_sw, T, initial_values[0], final_values[0],
            initial_values[1], final_values[1], initial_values[2], final_values[2],
            final_values[3], final_values[4], final_values[5]);

    file.close();
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
