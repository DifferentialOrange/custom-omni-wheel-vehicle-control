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

    QVector<double> t_noeps;
    QVector<double> nu_1_noeps;
    QVector<double> nu_2_noeps;
    QVector<double> nu_3_noeps;
    QVector<double> x_noeps;
    QVector<double> y_noeps;
    QVector<double> theta_noeps;
    QVector<double> chi_2_noeps;

    QVector<double> t;
    QVector<double> nu_1;
    QVector<double> nu_2;
    QVector<double> nu_3;
    QVector<double> x;
    QVector<double> y;
    QVector<double> theta;
    QVector<double> chi_2;

    double T;

    bool plotted;
    QCPCurve *trajectory_noeps;
    QCPCurve *trajectory_eps;
};

#endif // MAINWINDOW_H
