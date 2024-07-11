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
