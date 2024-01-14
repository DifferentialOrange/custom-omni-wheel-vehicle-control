#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "vector.h"
#include "qcustomplot.h"

namespace Ui
{
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();



private slots:
    void on_pushButton_compute_clicked();

private:
    Ui::MainWindow *ui;

    QVector<double> t_symm;
    QVector<double> nu_1_symm;
    QVector<double> nu_2_symm;
    QVector<double> nu_3_symm;
    QVector<double> x_symm;
    QVector<double> y_symm;
    QVector<double> theta_symm;

    QVector<double> t;
    QVector<double> nu_1;
    QVector<double> nu_2;
    QVector<double> nu_3;
    QVector<double> x;
    QVector<double> y;
    QVector<double> theta;

    Vector<6> initial_values, final_values;
    double t_sw, T;

    QVector<double> P_real, P_advice;
    QVector<double> N_1, N_2, N_3;

    bool plotted;
    QCPCurve *trajectory_minus_symm;
    QCPCurve *trajectory_plus_symm;
    QCPCurve *trajectory_minus;
    QCPCurve *trajectory_plus;
};

#endif // MAINWINDOW_H
