#include "parameters.h"
#include "vector.h"
#include <QVector>
#include <QDebug>
#include <gsl/gsl_linalg.h>

// dot x, dot y, dot theta, dot chi 1, dot chi 2, dot chi 3, x, y, theta, chi 1, chi 2, chi 3

//double f_margin = 1e-4;

std::array<double, 3> compute_N(double d_chi_1, double d_chi_2, double d_chi_3,
                                double d_theta, double U_1, double U_2, double U_3,
                                double sign_v_n_1, double sign_v_n_2, double sign_v_n_3)
{
    double g = 9.81;
    // works only with certain conditions:
    // S = Q (Delta = 0), corrent for symmetrical vehicle
    double a_data[] = {1, 1, 1,
                       parameters::symmetrical::Delta - parameters::symmetrical::delta[0] * sin(parameters::symmetrical::alpha[0])
                       - parameters::mu_n * sign_v_n_1 * sin(parameters::symmetrical::beta[0]),
                       parameters::symmetrical::Delta - parameters::symmetrical::delta[1] * sin(parameters::symmetrical::alpha[1])
                       - parameters::mu_n * sign_v_n_2 * sin(parameters::symmetrical::beta[1]),
                       parameters::symmetrical::Delta - parameters::symmetrical::delta[2] * sin(parameters::symmetrical::alpha[2])
                       - parameters::mu_n * sign_v_n_3 * sin(parameters::symmetrical::beta[2]),
                       parameters::symmetrical::delta[0] * cos(parameters::symmetrical::alpha[0])
                       + parameters::mu_n * sign_v_n_1 * cos(parameters::symmetrical::beta[0]),
                       parameters::symmetrical::delta[1] * cos(parameters::symmetrical::alpha[1])
                       + parameters::mu_n * sign_v_n_2 * cos(parameters::symmetrical::beta[1]),
                       parameters::symmetrical::delta[2] * cos(parameters::symmetrical::alpha[2])
                       + parameters::mu_n * sign_v_n_3 * cos(parameters::symmetrical::beta[2])};
    double b_data[] = {g,
                       parameters::lambda * parameters::lambda * (
                            d_chi_1 * d_theta * sin(parameters::symmetrical::beta[0]) +
                            d_chi_2 * d_theta * sin(parameters::symmetrical::beta[1]) +
                            d_chi_3 * d_theta * sin(parameters::symmetrical::beta[2])
                       ) -  cos(parameters::symmetrical::beta[0]) * (parameters::c1 * U_1 - parameters::c2 * d_chi_1) -
                            cos(parameters::symmetrical::beta[1]) * (parameters::c1 * U_2 - parameters::c2 * d_chi_2) -
                            cos(parameters::symmetrical::beta[2]) * (parameters::c1 * U_3 - parameters::c2 * d_chi_3),
                       parameters::lambda * parameters::lambda * (
                            - d_chi_1 * d_theta * cos(parameters::symmetrical::beta[0]) -
                            d_chi_2 * d_theta * cos(parameters::symmetrical::beta[1]) -
                            d_chi_3 * d_theta * cos(parameters::symmetrical::beta[2])
                       ) -  sin(parameters::symmetrical::beta[0]) * (parameters::c1 * U_1 - parameters::c2 * d_chi_1) -
                       sin(parameters::symmetrical::beta[1]) * (parameters::c1 * U_2 - parameters::c2 * d_chi_2) -
                       sin(parameters::symmetrical::beta[2]) * (parameters::c1 * U_3 - parameters::c2 * d_chi_3)
                      };


    gsl_matrix_view A = gsl_matrix_view_array(a_data, 3, 3);
    gsl_vector_view b = gsl_vector_view_array(b_data, 3);
    gsl_vector*     x = gsl_vector_alloc(3);

    gsl_permutation* p = gsl_permutation_alloc(3);
    int s;
    gsl_linalg_LU_decomp(&A.matrix, p, &s);
    gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);

    return { x->data[0], x->data[1], x->data[2] };
}

double f_sign(double a)
{
    return atan(1e4 * a) * 2 / M_PI;
//    if (fabs(a) < f_margin)
//        return a / f_margin;
//    else if (a < 0)
//        return -1;
//    else
//        return 1;
}

double v_m_i_n_sign(int i, double dot_x, double dot_y, double dot_theta, double dot_chi, double theta)
{
    return f_sign(
        dot_x * cos(theta + parameters::symmetrical::beta[i]) +
        dot_y * sin(theta + parameters::symmetrical::beta[i]) +
        dot_theta * (
            parameters::symmetrical::delta[i] *
            sin(parameters::symmetrical::beta[i] - parameters::symmetrical::alpha[i]) +
            parameters::symmetrical::Delta * cos(parameters::symmetrical::beta[i])
        )
    );
}

double v_m_i_tau_sign(int i, double dot_x, double dot_y, double dot_theta, double dot_chi, double theta)
{
    return f_sign(
        - dot_x * sin(theta + parameters::symmetrical::beta[i]) +
        dot_y * cos(theta + parameters::symmetrical::beta[i]) +
        dot_theta * (
            parameters::symmetrical::delta[i] *
            cos(parameters::symmetrical::beta[i] - parameters::symmetrical::alpha[i]) -
            parameters::symmetrical::Delta * sin(parameters::symmetrical::beta[i])
        ) - dot_chi
    );
}

Vector<12> rightpart(double t, Vector<12> x, Vector<3> control_minus, Vector<3> control_plus, double t_sw)
{
    Vector<12> rez;
    Vector<3> u;

    if (t < t_sw)
        u = control_minus;
    else
        u = control_plus;

    std::array<double, 3> v_m_i_n = {
        v_m_i_n_sign(0, x[0], x[1], x[2], x[3], x[8]),
        v_m_i_n_sign(1, x[0], x[1], x[2], x[4], x[8]),
        v_m_i_n_sign(2, x[0], x[1], x[2], x[5], x[8]),
    };

    std::array<double, 3> v_m_i_tau = {
        v_m_i_tau_sign(0, x[0], x[1], x[2], x[3], x[8]),
        v_m_i_tau_sign(1, x[0], x[1], x[2], x[4], x[8]),
        v_m_i_tau_sign(2, x[0], x[1], x[2], x[5], x[8]),
    };

    static std::array<double, 3> N = compute_N(x[3], x[4], x[5], x[2], u[0], u[1], u[2], v_m_i_tau[0], v_m_i_tau[1], v_m_i_tau[2]);

    rez[0] =
        - parameters::mu_n * (
            v_m_i_n[0] * cos(x[8] + parameters::symmetrical::beta[0]) * N[0] +
            v_m_i_n[1] * cos(x[8] + parameters::symmetrical::beta[1]) * N[1] +
            v_m_i_n[2] * cos(x[8] + parameters::symmetrical::beta[2]) * N[2]
        ) + parameters::mu_tau * (
            v_m_i_tau[0] * sin(x[8] + parameters::symmetrical::beta[0]) * N[0] +
            v_m_i_tau[1] * sin(x[8] + parameters::symmetrical::beta[1]) * N[1] +
            v_m_i_tau[2] * sin(x[8] + parameters::symmetrical::beta[2]) * N[2]
        )
    ;

    rez[1] =
        - parameters::mu_n * (
            v_m_i_n[0] * sin(x[8] + parameters::symmetrical::beta[0]) * N[0] +
            v_m_i_n[1] * sin(x[8] + parameters::symmetrical::beta[1]) * N[1] +
            v_m_i_n[2] * sin(x[8] + parameters::symmetrical::beta[2]) * N[2]
        ) - parameters::mu_tau * (
            v_m_i_tau[0] * cos(x[8] + parameters::symmetrical::beta[0]) * N[0] +
            v_m_i_tau[1] * cos(x[8] + parameters::symmetrical::beta[1]) * N[1] +
            v_m_i_tau[2] * cos(x[8] + parameters::symmetrical::beta[2]) * N[2]
        )
    ;

    rez[2] = (
        - parameters::mu_tau * (
                v_m_i_tau[0] * (
                    parameters::symmetrical::delta[0] *
                    cos(parameters::symmetrical::beta[0] - parameters::symmetrical::alpha[0]) -
                    parameters::symmetrical::Delta * sin(parameters::symmetrical::beta[0])
                ) * N[0] +
                v_m_i_tau[1] * (
                    parameters::symmetrical::delta[1] *
                    cos(parameters::symmetrical::beta[1] - parameters::symmetrical::alpha[1]) -
                    parameters::symmetrical::Delta * sin(parameters::symmetrical::beta[1])
                ) * N[1] +
                v_m_i_tau[2] * (
                    parameters::symmetrical::delta[2] *
                    cos(parameters::symmetrical::beta[2] - parameters::symmetrical::alpha[2]) -
                    parameters::symmetrical::Delta * sin(parameters::symmetrical::beta[2])
                ) * N[2]
        ) - parameters::mu_n * (
                v_m_i_n[0] * (
                    parameters::symmetrical::delta[0] *
                    sin(parameters::symmetrical::beta[0] - parameters::symmetrical::alpha[0]) +
                    parameters::symmetrical::Delta * cos(parameters::symmetrical::beta[0])
                ) * N[0] +
                v_m_i_n[1] * (
                    parameters::symmetrical::delta[1] *
                    sin(parameters::symmetrical::beta[1] - parameters::symmetrical::alpha[1]) +
                    parameters::symmetrical::Delta * cos(parameters::symmetrical::beta[1])
                ) * N[1] +
                v_m_i_n[2] * (
                    parameters::symmetrical::delta[2] *
                    sin(parameters::symmetrical::beta[2] - parameters::symmetrical::alpha[2]) +
                    parameters::symmetrical::Delta * cos(parameters::symmetrical::beta[2])
                ) * N[2]
        )
    ) / parameters::symmetrical::Lambda /  parameters::symmetrical::Lambda;

    rez[3] = (
        parameters::mu_tau * v_m_i_tau[0] * N[0] + parameters::c1 * u[0] - parameters::c2 * x[3]
    ) / parameters::lambda /  parameters::lambda;

    rez[4] = (
        parameters::mu_tau * v_m_i_tau[1] * N[1] + parameters::c1 * u[1] - parameters::c2 * x[4]
    ) / parameters::lambda /  parameters::lambda;

    rez[5] = (
        parameters::mu_tau * v_m_i_tau[2] * N[2] + parameters::c1 * u[2] - parameters::c2 * x[5]
    ) / parameters::lambda /  parameters::lambda;

    rez[6] = x[0];
    rez[7] = x[1];
    rez[8] = x[2];
    rez[9] = x[3];
    rez[10] = x[4];
    rez[11] = x[5];

    return rez;
}

Vector<12> DOPRI8_friction_plot(double t_left, double t_right, Vector<6> initial_values,
                                Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                                QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                                QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                                QVector<double> &theta_vec, QVector<double> &d_chi_1, QVector<double> &d_chi_2, QVector<double> &d_chi_3, QVector<double> &d_theta,
                                QVector<double> &v_sign_tau_1, QVector<double> &v_sign_tau_2, QVector<double> &v_sign_tau_3,
                                QVector<double> &v_sign_n_1, QVector<double> &v_sign_n_2, QVector<double> &v_sign_n_3,
                                QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3)
{
    double h = (t_right - t_left) / 1e7;
    double h_new;
    bool switch_flag = false;
    bool last_flag = false;

    // dot x, dot y, dot theta, dot chi 1, dot chi 2, dot chi 3, x, y, theta, chi 1, chi 2, chi 3
    double tl = t_left;
    Vector<12> xl;
    xl[0] = initial_values[0] * cos(initial_values[5]) - initial_values[1] * sin(initial_values[5]);
    xl[1] = initial_values[0] * sin(initial_values[5]) + initial_values[1] * cos(initial_values[5]);
    xl[2] = initial_values[2] / parameters::symmetrical::Lambda;
    xl[3] = 0;
    xl[4] = 0;
    xl[5] = 0;
    xl[6] = initial_values[3];
    xl[7] = initial_values[4];
    xl[8] = initial_values[5];
    xl[9] = 0;
    xl[10] = 0;
    xl[11] = 0;

    Vector<12> stepx;

    Vector<12> errx;
    double err, coef;
    double coefmax = 5;

    Vector<12> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

    t_vec.append(tl);
//    nu1_vec.append(xl[0]);
//    nu2_vec.append(xl[1]);
//    nu3_vec.append(xl[2]);
    x_vec.append(xl[6]);
    y_vec.append(xl[7]);
    theta_vec.append(xl[8]);
    d_chi_1.append(0);
    d_chi_2.append(0);
    d_chi_3.append(0);
    d_theta.append(0);

    std::array<double, 3> v_m_i_tau = {
        v_m_i_tau_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]),
        v_m_i_tau_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]),
        v_m_i_tau_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]),
    };
    auto N = compute_N(xl[3], xl[4], xl[5], xl[2], control_minus[0], control_minus[1], control_minus[2], v_m_i_tau[0], v_m_i_tau[1], v_m_i_tau[2]);
    N_1.append(N[0]);
    N_2.append(N[1]);
    N_3.append(N[2]);

//    v_sign_n_1.append(v_m_i_n_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]));
//    v_sign_n_2.append(v_m_i_n_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]));
//    v_sign_n_3.append(v_m_i_n_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]));

    v_sign_n_1.append(xl[0]);
    v_sign_n_2.append(xl[1]);
    v_sign_n_3.append(xl[2]);

    v_sign_tau_1.append(v_m_i_tau_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]));
    v_sign_tau_2.append(v_m_i_tau_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]));
    v_sign_tau_3.append(v_m_i_tau_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]));

    int i = 0;

    while (tl + h < t_right || last_flag)
    {

        //switch point
        if (tl < t_sw && tl + h > t_sw && !switch_flag)
        {
            h = t_sw - tl;
            switch_flag = true;
        }

        k1 = rightpart(tl, xl, control_minus, control_plus, t_sw);

        k2 = rightpart(tl + h / 18, xl + h * k1 / 18, control_minus, control_plus, t_sw);

        k3 = rightpart(tl + h / 12, xl + h * (k1  / 48 + k2 / 16) , control_minus, control_plus, t_sw);

        k4 = rightpart(tl + h / 8, xl + h * (k1  / 32 + k3 * 3 / 32) , control_minus, control_plus, t_sw);

        k5 = rightpart(tl + h * 5 / 16, xl + h * (k1  * 5 / 16 + k3 * (-75) / 64 + k4 * 75 / 64) , control_minus, control_plus, t_sw);

        k6 = rightpart(tl + h * 3 / 8, xl + h * (k1  * 3 / 80 + k4 * 3 / 16 + k5 * 3 / 20) , control_minus, control_plus, t_sw);

        k7 = rightpart(tl + h * 59. / 400., xl + h * (k1  * (29443841. / 614563906.) + k4 * (77736538. / 692538347.) + k5 * ((-28693883.) / 1125000000.) + 
                k6 * (23124283. / 1800000000.) ), control_minus, control_plus, t_sw);

        k8 = rightpart(tl + h * 93. / 200., xl + h * (k1  * (16016141. / 946692911.) + k4 * (61564180. / 158732637.) + k5 * (22789713. / 633445777.) + 
                k6 * (545815736. / 2771057229.) + k7 * ((-180193667.) / 1043307555.) ), control_minus, control_plus, t_sw);

        k9 = rightpart(tl + h * (5490023248. / 9719169821.), xl + h * (k1  * (39632708. / 573591083.) + k4 * ((-433636366.) / 683701615.) + k5 * ((-421739975.) / 2616292301.) +
                k6 * (100302831. / 723423059.) + k7 * (790204164. / 839813087.) + k8 * (800635310. / 3783071287.) ), control_minus, control_plus, t_sw);

        k10 = rightpart(tl + h * 13 / 20, xl + h * (k1  * (246121993. / 1340847787.) + k4 * ((-37695042795.) / 15268766246.) + k5 * ((-309121744.) / 1061227803.) +
                k6 * ((-12992083.) / 490766935.) + k7 * (6005943493. / 2108947869.) + k8 * (393006217. / 1396673457) + k9 * (123872331. / 1001029789.) ), control_minus, control_plus, t_sw);

        k11 = rightpart(tl + h * (1201146811. / 1299019798.), xl + h * (k1  * ((-1028468189.) / 846180014.) + k4 * (8478235783. / 508512852.) + k5 * (1311729495. / 1432422823.) +
                k6 * ((-10304129995.) / 1701304382.) + k7 * ((-48777925059.) / 3047939560.) + k8 * (15336726248. / 1032824649.) + k9 * ((-45442868181.) / 3398467696.) + 
                k10 * (3065993473. / 597172653.) ), control_minus, control_plus, t_sw);

        k12 = rightpart(tl + h, xl + h * (k1  * (185892177. / 718116043.) + k4 * ((-3185094517.) / 667107341.) + k5 * ((-477755414.) / 1098053517.) + 
                 k6 * ((-703635378.) / 230739211.) + k7 * (5731566787. / 1027545527.) + k8 * (5232866602. / 850066563.) + k9 * ((-4093664535.) / 808688257.) + 
                k10 * (3962137247. / 1805957418.) + k11 * (65686358. / 487910083.) ), control_minus, control_plus, t_sw);


        k13 = rightpart(tl + h, xl + h * (k1  * (403863854. / 491063109.) + k4 * ((-5068492393.) / 434740067.) + k5 * ((-411421997.) / 543043805.) + 
                 k6 * (652783627. / 914296604.) + k7 * (11173962825. / 925320556.) + k8 * ((-13158990841.) / 6184727034.) + k9 * (3936647629. / 1978049680.) + 
                k10 * ((-160528059.) / 685178525.) + k11 * (248638103. / 1413531060.) ), control_minus, control_plus, t_sw);


        stepx = h * (k1 * (14005451. / 335480064.) + k6 * ((-59238493.) / 1068277825.) + k7 * (181606767. / 758867731.) + k8 * (561292985. / 797845732.) + 
                k9 * ((-1041891430.) / 1371343529.) + k10 * (760417239. / 1151165299.) + k11 * (118820643. / 751138087.) + k12 * ((-528747749.) / 2220607170.) + k13 / 4 );


        errx = h * (k1 * (13451932. / 455176623.) + k6 * ((-808719846.) / 976000145.) + k7 * (1757004468. / 5645159321.) + k8 * (656045339. / 265891186.) + 
                k9 * ((-3867574721.) / 1518517206.) + k10 * (465885868. / 322736535.) + k11 * (53011238. / 667516719.) + k12 * 2 / 45 );


        err = (stepx - errx).euclidDistance();

        //step correction section
        if (err < 1e-15)
        {
            tl += h;
            xl += stepx;

            t_vec.append(tl);
//            nu1_vec.append(xl[0]);
//            nu2_vec.append(xl[1]);
//            nu3_vec.append(xl[2]);
            x_vec.append(xl[6]);
            y_vec.append(xl[7]);
            theta_vec.append(xl[8]);

            d_chi_1.append(xl[3]);
            d_chi_2.append(xl[4]);
            d_chi_3.append(xl[5]);
            d_theta.append(xl[2]);
            std::array<double, 3> v_m_i_tau = {
                v_m_i_tau_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]),
                v_m_i_tau_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]),
                v_m_i_tau_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]),
            };
            auto N = compute_N(xl[3], xl[4], xl[5], xl[2], control_minus[0], control_minus[1], control_minus[2], v_m_i_tau[0], v_m_i_tau[1], v_m_i_tau[2]);
            N_1.append(N[0]);
            N_2.append(N[1]);
            N_3.append(N[2]);

            v_sign_n_1.append(v_m_i_n_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]));
            v_sign_n_2.append(v_m_i_n_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]));
            v_sign_n_3.append(v_m_i_n_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]));

            v_sign_tau_1.append(v_m_i_tau_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]));
            v_sign_tau_2.append(v_m_i_tau_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]));
            v_sign_tau_3.append(v_m_i_tau_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]));

            coefmax = 5;

            h *= 2; 
        }
        else
        {
            coef = pow(precision::DOPRI8_error_EPS / err, 1.0 / (8 + 1));

            if (coef > coefmax)
                coef = coefmax;
            else if (coef < 0.2)
                coef = 0.2;

            h_new = 0.9 * h * coef;


            if (precision::DOPRI8_error_EPS > err)
            {           
                tl += h;
                xl += stepx;

                t_vec.append(tl);
//                nu1_vec.append(xl[0]);
//                nu2_vec.append(xl[1]);
//                nu3_vec.append(xl[2]);
                x_vec.append(xl[6]);
                y_vec.append(xl[7]);
                theta_vec.append(xl[8]);

                d_chi_1.append(xl[3]);
                d_chi_2.append(xl[4]);
                d_chi_3.append(xl[5]);
                d_theta.append(xl[2]);
                std::array<double, 3> v_m_i_tau = {
                    v_m_i_tau_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]),
                    v_m_i_tau_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]),
                    v_m_i_tau_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]),
                };
                auto N = compute_N(xl[3], xl[4], xl[5], xl[2], control_minus[0], control_minus[1], control_minus[2], v_m_i_tau[0], v_m_i_tau[1], v_m_i_tau[2]);
                N_1.append(N[0]);
                N_2.append(N[1]);
                N_3.append(N[2]);

                v_sign_n_1.append(v_m_i_n_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]));
                v_sign_n_2.append(v_m_i_n_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]));
                v_sign_n_3.append(v_m_i_n_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]));

                v_sign_tau_1.append(v_m_i_tau_sign(0, xl[0], xl[1], xl[2], xl[3], xl[8]));
                v_sign_tau_2.append(v_m_i_tau_sign(1, xl[0], xl[1], xl[2], xl[4], xl[8]));
                v_sign_tau_3.append(v_m_i_tau_sign(2, xl[0], xl[1], xl[2], xl[5], xl[8]));

                coefmax = 5;
            }
            else
                coefmax = 1;

            h = h_new;
        }

        //for correct plotting; there could be too few points
        if (h > (t_right - t_left) / 100)
        {
            h = (t_right - t_left) / 100;
        }

        //last step correction
        if (tl + h >= t_right && last_flag == false)
        {
            h = t_right - tl;
            last_flag = true;
        }
        else
            last_flag = false;

        i++;

        if (i % 1000 == 0) {
            qDebug() << "iteration " << i << '\n';
            qDebug() << "tl " << tl << '\n';
            qDebug() << "h " << h << '\n';
            qDebug() << "err " << err << '\n';
            qDebug() << "xl " << xl[0] << ' ' <<  xl[1] << ' ' << xl[2] << '\n';
        }
    }

    return xl;
}
