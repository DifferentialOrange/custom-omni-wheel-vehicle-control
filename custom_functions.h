#pragma once

#include <cmath>
#include <array>
#include "parameters.h"
#include "vector.h"

std::array<double, 9> matr_sum (std::array<double, 9> A, std::array<double, 9> B);
std::array<double, 9> matr_mult (std::array<double, 9> A, std::array<double, 9> B);
std::array<double, 9> scal_mult (std::array<double, 9> A, double l);
Vector<3> vec_mult (std::array<double, 9> A, Vector<3> v);
std::array<double, 9> matr_T (std::array<double, 9> A);
double matr_det (std::array<double, 9> A);
std::array<double, 9> matr_rev (std::array<double, 9> A);
static const std::array<double, 9> matr_E({1, 0, 0, 0, 1, 0, 0, 0, 1});
std::array<double, 9> Sigma_matr (std::array<double, 3> alpha, std::array<double, 3> beta,
                                  std::array<double, 3> delta, double Delta, double Lambda);
std::array<double, 9> A_matr (std::array<double, 3> alpha, std::array<double, 3> beta,
                              std::array<double, 3> delta, double Delta, double Lambda);
