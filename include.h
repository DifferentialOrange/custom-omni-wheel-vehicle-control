#pragma once

#include <cmath>
#include <functional>
#include <cstdio>
#include <iostream>
#include <random>
#include <algorithm>
#include <QVector>

#include <gsl/gsl_linalg.h>

#include "vector.h"
#include "parameters.h"

Vector<6> DOPRI8_symmetrical_plot(double t_left, double t_right, Vector<6> initial_values,
                    Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec,
                    QVector<double> &P_real, QVector<double> &P_advice, QVector<double> &PT_advice,
                    QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3,
                    QVector<double> &U1, QVector<double> &U2, QVector<double> &U3);
