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

void DOPRI8_symmetrical_plot(double t_left, double t_right,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec, QVector<double> &chi_2_vec, double eps);
