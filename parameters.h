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

namespace parameters
{
    const double lambda  = sqrt(5e-4);

    const double c1      = 1e-2;
    const double c2      = 2.5e-4;

    namespace symmetrical
    {
        const double Lambda = sqrt(5e-2);

        const double m = 3.0;
        const double rho     = 0.1;
        const double r       = 0.05;

        const double A1      = m + 3 * lambda * lambda / 2 / r / r;
        const double A2      = m + 3 * lambda * lambda / 2 / r / r;
        const double A3      = 1 + 3 * rho * rho * lambda * lambda / (Lambda * Lambda * r * r);
        const double L       = m / (Lambda * A1); //L1 = L2 == L

        const double kappa   = 3 * c2 / (2 * A1 * r * r); //kappa1 == kappa2 == kappa
        const double kappa3  = 3 * c2 * rho * rho / (Lambda * Lambda * A3 * r * r);

        const std::array<double, 3> alpha = {- M_PI / 6, M_PI / 2, 7 * M_PI / 6};
        const std::array<double, 3> beta = alpha;
        const std::array<double, 3> delta = {rho, rho, rho};
        const double Delta = 0;
    }

    namespace final
    {
        const double Lambda = 5.0 / 4;
        const double Delta = 0.5;

        const std::array<double, 3> alpha = {0, M_PI / 2, M_PI};
        const std::array<double, 3> beta = alpha;
        const std::array<double, 3> delta = {5.0 / 4, 2.0, 5.0 / 4};
    }
}
