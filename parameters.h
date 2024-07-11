#pragma once
#include <array>
#include <cmath>

namespace precision
{
    const double DOPRI8_error_EPS = 1e-9;
    const double Newton_EPS = 1e-9 ;
    const double double_EPS = 1e-13;
    const double hder = 1e-4; //numerical derivative step

    const double angle_step = M_PI / 12;
    const double N_step = 0.1;
}

namespace parameters2
{
    const double m = 3;
    const double R = 0.05;
    const double h = 0.2;
    const double rho = 0.15;
    const double Lambda = sqrt(0.05);
    const double lambda = sqrt(5e-4);
    const double c_1 = 0.01;
    const double c_2 = 2.5 * 1e-4;
    const std::array<double, 3> alpha = {- M_PI / 6, M_PI / 2, 7 * M_PI / 6};
    const std::array<double, 3> beta = alpha;
    const double g = 9.81;
}
