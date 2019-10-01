#include "include.h"
#include "vector.h"
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <QDebug>

std::pair<bool, Vector<6>> control_solve(double t_left, double t_right, Vector<6> initial_values, Vector<6> final_values,
                        Vector<6> control, double t_sw,
                        std::array<double, 3> alpha, std::array<double, 3> beta,
                        std::array<double, 3> delta, double Delta, double Lambda)
{
    auto current_values = DOPRI8_custom(t_left, t_right, initial_values, {control[0], control[1], control[2]},
                                {control[3], control[4], control[5]}, t_sw, alpha, beta, delta, Delta, Lambda);

    auto start_dev = (final_values - current_values).euclidDistance();

    while ((final_values - current_values).euclidDistance() >
           precision::Newton_EPS * std::max(1.0, final_values.euclidDistance()))
    {
        double der_raw[36];
        for (unsigned i = 0; i < 6; ++i)
        {
            Vector<6> control_step = {0, 0, 0, 0, 0, 0};
            control_step[i] = precision::hder;

            auto plus_step = control + control_step;
            auto plus_sol = DOPRI8_custom(t_left, t_right, initial_values, {plus_step[0], plus_step[1], plus_step[2]},
            {plus_step[3], plus_step[4], plus_step[5]}, t_sw, alpha, beta, delta, Delta, Lambda);
            auto minus_step = control - control_step;
            auto minus_sol = DOPRI8_custom(t_left, t_right, initial_values, {minus_step[0], minus_step[1], minus_step[2]},
            {minus_step[3], minus_step[4], minus_step[5]}, t_sw, alpha, beta, delta, Delta, Lambda);

            auto line = (plus_sol - minus_sol) / 2.0 / precision::hder;

            for (unsigned j = 0; j < 6; ++j)
                der_raw[i + j * 6] = line[j];
        }

        gsl_matrix_view der_matr = gsl_matrix_view_array(der_raw, 6, 6);

        gsl_permutation* p = gsl_permutation_alloc(6);
        int s;
        gsl_linalg_LU_decomp(&der_matr.matrix, p, &s);

        double inv_raw[36];
        gsl_matrix_view der_inv = gsl_matrix_view_array(inv_raw, 6, 6);

        gsl_linalg_LU_invert (&der_matr.matrix, p, &der_inv.matrix);

        auto f_V = current_values - final_values;
        double f_raw[] = {f_V[0], f_V[1], f_V[2],
                         f_V[3], f_V[4], f_V[5]};

        gsl_matrix_view f_vec = gsl_matrix_view_array(f_raw, 6, 1);

        double mult_raw[6];
        gsl_matrix_view mult_vec = gsl_matrix_view_array(mult_raw, 6, 1);

        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, &der_inv.matrix, &f_vec.matrix,
                        0.0, &mult_vec.matrix);

        for (unsigned j = 0; j < 6; ++j)
            control[j] -= mult_vec.matrix.data[j];

        current_values = DOPRI8_custom(t_left, t_right, initial_values, {control[0], control[1], control[2]},
                                    {control[3], control[4], control[5]}, t_sw, alpha, beta, delta, Delta, Lambda);

        if ((final_values - current_values).euclidDistance() > start_dev * 10)
        {
            qDebug() << "Bad convergence\n";

            return {false, control};
        }
        qDebug() << (final_values - current_values).euclidDistance() << '\n';
    }

    return {true, control};
}
