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
#include "symm_functions.h"

double DOPRI8_par_integrate(double t_left, double t_right,
                            std::function<double(double, double, double, double, double)>,
                            double D3, double W3, double D6);

Vector<6> predict_control(double t_sw, double T, double nu_1_0, double nu_1_T,
                          double nu_2_0, double nu_2_T, double nu_3_0, double nu_3_T,
                          double x_T, double y_T, double theta_0, double theta_T);

std::pair<bool, Vector<6>> control_solve(double t_left, double t_right, Vector<6> initial_values,
                        Vector<6> final_values, Vector<6> control, double t_sw,
                        std::array<double, 3> alpha, std::array<double, 3> beta,
                        std::array<double, 3> delta, double Delta, double Lambda);

Vector<6> custom_control_find(Vector<6> control, double t_sw, double T, Vector<6> initial_values,
                              Vector<6> final_values);

Vector<6> DOPRI8_symmetrical_plot(double t_left, double t_right, Vector<6> initial_values,
                    Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec,
                    QVector<double> &P_real, QVector<double> &P_advice,
                    QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3,
                    QVector<double> &U1, QVector<double> &U2, QVector<double> &U3);

Vector<6> DOPRI8_custom(double t_left, double t_right, Vector<6> initial_values,
                    Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                    std::array<double, 3> alpha, std::array<double, 3> beta,
                    std::array<double, 3> delta, double Delta, double Lambda);

Vector<6> DOPRI8_final_plot(double t_left, double t_right, Vector<6> initial_values,
                    Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec);
