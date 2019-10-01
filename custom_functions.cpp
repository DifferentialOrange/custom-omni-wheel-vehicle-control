#include <cmath>
#include <array>
#include "parameters.h"
#include "vector.h"
#include "custom_functions.h"

std::array<double, 9> matr_sum (std::array<double, 9> A, std::array<double, 9> B)
{
    return {A[0] + B[0], A[1] + B[1], A[2] + B[2],
            A[3] + B[3], A[4] + B[4], A[5] + B[5],
            A[6] + B[6], A[7] + B[7], A[8] + B[8]};
}

std::array<double, 9> matr_mult (std::array<double, 9> A, std::array<double, 9> B)
{
    return {A[0] * B[0] + A[1] * B[3] + A[2] * B[6], A[0] * B[1] + A[1] * B[4] + A[2] * B[7],
            A[0] * B[2] + A[1] * B[5] + A[2] * B[8], A[3] * B[0] + A[4] * B[3] + A[5] * B[6],
            A[3] * B[1] + A[4] * B[4] + A[5] * B[7], A[3] * B[2] + A[4] * B[5] + A[5] * B[8],
            A[6] * B[0] + A[7] * B[3] + A[8] * B[6], A[6] * B[1] + A[7] * B[4] + A[8] * B[7],
            A[6] * B[2] + A[7] * B[5] + A[8] * B[8]};
}

std::array<double, 9> scal_mult (std::array<double, 9> A, double l)
{
    return {A[0] * l, A[1] * l, A[2] * l,
            A[3] * l, A[4] * l, A[5] * l,
            A[6] * l, A[7] * l, A[8] * l};
}

Vector<3> vec_mult (std::array<double, 9> A, Vector<3> v)
{
    return {A[0] * v[0] + A[1] * v[1] + A[2] * v[2],
            A[3] * v[0] + A[4] * v[1] + A[5] * v[2],
            A[6] * v[0] + A[7] * v[1] + A[8] * v[2]};
}

std::array<double, 9> matr_T (std::array<double, 9> A)
{
    return {A[0], A[3], A[6],
            A[1], A[4], A[7],
            A[2], A[5], A[8]};
}

double matr_det (std::array<double, 9> A)
{
    return A[0] * A[4] * A[8] + A[2] * A[3] * A[7] + A[1] * A[5] * A[6] -
            A[2] * A[4] * A[6] - A[1] * A[3] * A[8] - A[0] * A[5] * A[7];
}
std::array<double, 9> matr_rev (std::array<double, 9> A)
{
    return {(A[4] * A[8] - A[5] * A[7]) / matr_det(A), (A[2] * A[7] - A[1] * A[8]) / matr_det(A),
            (A[1] * A[5] - A[2] * A[4]) / matr_det(A), (A[5] * A[6] - A[3] * A[8]) / matr_det(A),
            (A[0] * A[8] - A[2] * A[6]) / matr_det(A), (A[2] * A[3] - A[0] * A[5]) / matr_det(A),
            (A[3] * A[7] - A[4] * A[6]) / matr_det(A), (A[1] * A[6] - A[0] * A[7]) / matr_det(A),
            (A[0] * A[4] - A[1] * A[3]) / matr_det(A)};
}

std::array<double, 9> Sigma_matr (std::array<double, 3> alpha, std::array<double, 3> beta,
                                  std::array<double, 3> delta, double Delta, double Lambda)
{
    return {-sin(beta[0]), cos(beta[0]), -Delta * sin(beta[0]) / Lambda
                                    + delta[0] * cos(beta[0] - alpha[0]) / Lambda,
            -sin(beta[1]), cos(beta[1]), -Delta * sin(beta[1]) / Lambda
                                    + delta[1] * cos(beta[1] - alpha[1]) / Lambda,
            -sin(beta[2]), cos(beta[2]), -Delta * sin(beta[2]) / Lambda
                                    + delta[2] * cos(beta[2] - alpha[2]) / Lambda};
}

std::array<double, 9> A_matr (std::array<double, 3> alpha, std::array<double, 3> beta,
                              std::array<double, 3> delta, double Delta, double Lambda)
{
    auto Sigma = Sigma_matr(alpha, beta, delta, Delta, Lambda);
    return matr_sum(matr_E, scal_mult(matr_mult(matr_T(Sigma), Sigma), parameters::lambda * parameters::lambda));
}
