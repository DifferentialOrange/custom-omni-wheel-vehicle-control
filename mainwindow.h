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
    void on_pushButton_generate_clicked();

private:
    Ui::MainWindow *ui;

    QVector<double> t;
    QVector<double> nu_1;
    QVector<double> nu_2;
    QVector<double> nu_3;
    QVector<double> x;
    QVector<double> y;
    QVector<double> theta;

    QVector<double> t_err;
    QVector<double> nu_1_err;
    QVector<double> nu_2_err;
    QVector<double> nu_3_err;
    QVector<double> x_err;
    QVector<double> y_err;
    QVector<double> theta_err;

    Vector<6> initial_values, final_values;
    double t_sw, T;

    bool plotted;

    QCPCurve *trajectory_minus;
    QCPCurve *trajectory_plus;
    QCPCurve *trajectory_minus_err;
    QCPCurve *trajectory_plus_err;
};

#endif // MAINWINDOW_H
