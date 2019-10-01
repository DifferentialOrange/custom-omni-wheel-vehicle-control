#include "include.h"
#include <gsl/gsl_linalg.h>

Vector<6> predict_control(double t_sw, double T, double nu_1_0, double nu_1_T, double nu_2_0, double nu_2_T, double nu_3_0, double nu_3_T, 
                            double x_T, double y_T, double theta_T)
{

    double a_first_data[] = {    1, 0, 1 / parameters::symmetrical::kappa3, 0, 0, 0,
                                exp(-parameters::symmetrical::kappa3 * t_sw), -exp(-parameters::symmetrical::kappa3 * t_sw), 1 / parameters::symmetrical::kappa3, -1 / parameters::symmetrical::kappa3, 0, 0,
                                0, exp(-parameters::symmetrical::kappa3 * T), 0, 1 / parameters::symmetrical::kappa3, 0, 0,
                                -1 / (parameters::symmetrical::kappa3 * parameters::symmetrical::Lambda), 0, 0, 0, 1, 0,
                                -exp(-parameters::symmetrical::kappa3 * t_sw) / (parameters::symmetrical::kappa3 * parameters::symmetrical::Lambda), exp(-parameters::symmetrical::kappa3 * t_sw) / (parameters::symmetrical::kappa3 * parameters::symmetrical::Lambda), t_sw / (parameters::symmetrical::kappa3 * parameters::symmetrical::Lambda), -t_sw / (parameters::symmetrical::kappa3 * parameters::symmetrical::Lambda), 1, -1,
                                0, -exp(-parameters::symmetrical::kappa3 * T) / (parameters::symmetrical::kappa3 * parameters::symmetrical::Lambda), 0, T / (parameters::symmetrical::kappa3 * parameters::symmetrical::Lambda), 0, 1};
    double b_first_data[] = {nu_3_0, 0, nu_3_T, 0, 0, theta_T};


    gsl_matrix_view A_first = gsl_matrix_view_array(a_first_data, 6, 6);
    gsl_vector_view b_first = gsl_vector_view_array(b_first_data, 6);
    gsl_vector*     x_first = gsl_vector_alloc(6);

    gsl_permutation* p_first = gsl_permutation_alloc(6);
    int s;
    gsl_linalg_LU_decomp(&A_first.matrix, p_first, &s);
    gsl_linalg_LU_solve(&A_first.matrix, p_first, &b_first.vector, x_first);


    double D_3_minus = x_first->data[0];
    double D_3_plus  = x_first->data[1];
    double W_3_minus = x_first->data[2];
    double W_3_plus  = x_first->data[3];
    double D_6_minus = x_first->data[4];
    double D_6_plus  = x_first->data[5];

    double a_second_data[] = {  -1, gamma_s_plus_c(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, 0, gamma_s_minus_c(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, 0, 0,
                                0, -gamma_s_minus_c(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, -1, gamma_s_plus_c(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, 0, 0,
                                exp(-parameters::symmetrical::kappa * T) * cos(beta(t_sw, T, D_3_plus, W_3_plus, D_6_plus)), 0, gamma_s_plus_c(t_sw, T, D_3_plus, W_3_plus, D_6_plus), exp(-parameters::symmetrical::kappa * T) * sin(beta(t_sw, T, D_3_plus, W_3_plus, D_6_plus)), 0, gamma_s_minus_c(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 0, 0,
                                -exp(-parameters::symmetrical::kappa * T) * sin(beta(t_sw, T, D_3_plus, W_3_plus, D_6_plus)), 0, -gamma_s_minus_c(t_sw, T, D_3_plus, W_3_plus, D_6_plus), exp(-parameters::symmetrical::kappa * T) * cos(beta(t_sw, T, D_3_plus, W_3_plus, D_6_plus)), 0, gamma_s_plus_c(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 0, 0,
                                0, Gamma_W_1(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, 0, -Gamma_W_2(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, -1, 0,
                                0, Gamma_W_2(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, 0, Gamma_W_1(0, t_sw, D_3_minus, W_3_minus, D_6_minus), 0, 0, -1,
                                Gamma_D_1(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 0, Gamma_W_1(t_sw, T, D_3_plus, W_3_plus, D_6_plus), Gamma_D_2(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 0, -Gamma_W_2(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 1, 0,
                                -Gamma_D_2(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 0, Gamma_W_2(t_sw, T, D_3_plus, W_3_plus, D_6_plus), Gamma_D_1(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 0, Gamma_W_1(t_sw, T, D_3_plus, W_3_plus, D_6_plus), 0, 1};

    double b_second_data[] = {  -exp(-parameters::symmetrical::kappa * t_sw) * (cos(beta(0, t_sw, D_3_minus, W_3_minus, D_6_minus)) * nu_1_0 + sin(beta(0, t_sw, D_3_minus, W_3_minus, D_6_minus)) * nu_2_0),
                                -exp(-parameters::symmetrical::kappa * t_sw) * (-sin(beta(0, t_sw, D_3_minus, W_3_minus, D_6_minus)) * nu_1_0 + cos(beta(0, t_sw, D_3_minus, W_3_minus, D_6_minus)) * nu_2_0),
                                nu_1_T,
                                nu_2_T,
                                -nu_1_0 * Gamma_D_1(0, t_sw, D_3_minus, W_3_minus, D_6_minus) - nu_2_0 * Gamma_D_2(0, t_sw, D_3_minus, W_3_minus, D_6_minus),
                                nu_1_0 * Gamma_D_2(0, t_sw, D_3_minus, W_3_minus, D_6_minus) - nu_2_0 * Gamma_D_1(0, t_sw, D_3_minus, W_3_minus, D_6_minus),
                                x_T,
                                y_T};
    gsl_matrix_view A_second     = gsl_matrix_view_array(a_second_data, 8, 8);
    gsl_vector_view b_second     = gsl_vector_view_array(b_second_data, 8);
    gsl_vector*    x_second      = gsl_vector_alloc(8);

    gsl_permutation* p_second = gsl_permutation_alloc(8);
    gsl_linalg_LU_decomp(&A_second.matrix, p_second, &s);
    gsl_linalg_LU_solve(&A_second.matrix, p_second, &b_second.vector, x_second);

    double W_1_minus    = x_second->data[1];
    double W_1_plus     = x_second->data[2];
    double W_2_minus    = x_second->data[4];
    double W_2_plus     = x_second->data[5];


    Vector<6> control = {u_control_1(W_1_minus, W_2_minus, W_3_minus), u_control_2(W_1_minus, W_2_minus, W_3_minus), u_control_3(W_1_minus, W_2_minus, W_3_minus), 
                        u_control_1(W_1_plus, W_2_plus, W_3_plus), u_control_2(W_1_plus, W_2_plus, W_3_plus), u_control_3(W_1_plus, W_2_plus, W_3_plus)};

    gsl_permutation_free(p_first);
    gsl_permutation_free(p_second);
    gsl_vector_free(x_first);
    gsl_vector_free(x_second);

    return control;
}
