#include "parameters.h"
#include "vector.h"
#include <QVector>
#include <QDebug>
#include <math.h>
#include <gsl/gsl_linalg.h>

//nu1, nu2, nu3, x, y, theta

Vector<6> rightpart(double t, Vector<6> x, Vector<3> u)
{
    Vector<6> rez;

    rez[0] = (x[1] * x[2] / parameters::symmetrical::Lambda + parameters::c1 * (u[0] / 2 - u[1] + u[2] / 2)
            - 3 * parameters::c2 * x[0] / 2) / parameters::symmetrical::A1;
    rez[1] = (- x[0] * x[2] / parameters::symmetrical::Lambda + parameters::c1 * sqrt(3) * (u[0]  - u[2]) / 2
            - 3 * parameters::c2 * x[1] / 2) / parameters::symmetrical::A2;
    rez[2] = (parameters::c1 * parameters::symmetrical::rho * (u[0] + u[1] + u[2]) / parameters::symmetrical::Lambda
            - 3 * parameters::c2 * parameters::symmetrical::rho * parameters::symmetrical::rho * x[2]
            / (parameters::symmetrical::Lambda * parameters::symmetrical::Lambda)) / parameters::symmetrical::A3;
    rez[3] = x[0] * cos(x[5]) - x[1] * sin(x[5]);
    rez[4] = x[0] * sin(x[5]) + x[1] * cos(x[5]);
    rez[5] = x[2] / parameters::symmetrical::Lambda;

    return rez;
}

void compute_P(double t, Vector<6> x, Vector<3> u, QVector<double> &P_real, QVector<double> &P_advice)
{
    P_real.append(sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));

    double g = 9.81;

    double h = (1 + 3.0 * parameters::lambda * parameters::lambda / 2) * (x[0] * x[0] + x[1] * x[1]) +
            (1 + 3.0 * parameters::symmetrical::rho * parameters::symmetrical::rho * parameters::lambda * parameters::lambda /
             parameters::symmetrical::Lambda / parameters::symmetrical::Lambda) * x[2] * x[2];

    double good_advice = g * parameters::symmetrical::rho / parameters::c1 / sqrt(6) -
                (3 + sqrt(3)) * (parameters::lambda * parameters::lambda * h / sqrt(parameters::symmetrical::Lambda * parameters::symmetrical::Lambda +
                                                                                    3 * parameters::symmetrical::rho * parameters::symmetrical::rho *
                                                                                    parameters::lambda * parameters::lambda) +
                                 parameters::c2 * sqrt(h)) / 2 / parameters::c1 / sqrt(2 + 3 * parameters::lambda * parameters::lambda);

    P_advice.append(good_advice);
}

void compute_N(double t, Vector<6> x, Vector<3> control,
               QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3)
{
    double g = 9.81;

        //qDebug() << "u symm " << u[0] << " " << u[1] << " " << u[2] << '\n';

    double d_1 = parameters::symmetrical::delta[0];
    double d_2 = parameters::symmetrical::delta[1];
    double d_3 = parameters::symmetrical::delta[2];

    double D = parameters::symmetrical::Delta;

    double a_1 = parameters::symmetrical::alpha[0];
    double a_2 = parameters::symmetrical::alpha[1];
    double a_3 = parameters::symmetrical::alpha[2];
    double b_1 = parameters::symmetrical::beta[0];
    double b_2 = parameters::symmetrical::beta[1];
    double b_3 = parameters::symmetrical::beta[2];

    double L = parameters::symmetrical::Lambda;
    double l = parameters::lambda;

//    double D_1 = D * cos(b_1) + d_1 * sin(b_1 - a_1);
//    double D_2 = D * cos(b_2) + d_2 * sin(b_2 - a_2);
//    double D_3 = D * cos(b_3) + d_3 * sin(b_3 - a_3);

    double s_11 = -sin(b_1);
    double s_21 = -sin(b_2);
    double s_31 = -sin(b_3);
    double s_12 = cos(b_1);
    double s_22 = cos(b_2);
    double s_32 = cos(b_3);
    double s_13 = (- D * sin(b_1) + d_1 * cos(a_1 - b_1)) / L;
    double s_23 = (- D * sin(b_2) + d_2 * cos(a_2 - b_2)) / L;
    double s_33 = (- D * sin(b_3) + d_3 * cos(a_3 - b_3)) / L;
    double d_chi_1 = s_11 * x[0] + s_12 * x[1] + s_13 * x[2];
    double d_chi_2 = s_21 * x[0] + s_22 * x[1] + s_23 * x[2];
    double d_chi_3 = s_31 * x[0] + s_32 * x[1] + s_33 * x[2];
    double d_theta = x[2] / L;

    double a_data[] = {1, 1, 1,
                       - parameters::symmetrical::Delta + parameters::symmetrical::delta[0] * sin(parameters::symmetrical::alpha[0]),
                       - parameters::symmetrical::Delta + parameters::symmetrical::delta[1] * sin(parameters::symmetrical::alpha[1]),
                       - parameters::symmetrical::Delta + parameters::symmetrical::delta[2] * sin(parameters::symmetrical::alpha[2]),
                       - parameters::symmetrical::delta[0] * cos(parameters::symmetrical::alpha[0]),
                       - parameters::symmetrical::delta[1] * cos(parameters::symmetrical::alpha[1]),
                       - parameters::symmetrical::delta[2] * cos(parameters::symmetrical::alpha[2])};
    double b_data[] = {g,
                       - parameters::lambda * parameters::lambda * (
                            d_chi_1 * d_theta * cos(parameters::symmetrical::beta[0]) +
                            d_chi_2 * d_theta * cos(parameters::symmetrical::beta[1]) +
                            d_chi_3 * d_theta * cos(parameters::symmetrical::beta[2])
                       ) - sin(parameters::symmetrical::beta[0]) * (parameters::c1 * control[0] - parameters::c2 * d_chi_1) -
                           sin(parameters::symmetrical::beta[1]) * (parameters::c1 * control[1] - parameters::c2 * d_chi_2) -
                           sin(parameters::symmetrical::beta[2]) * (parameters::c1 * control[2] - parameters::c2 * d_chi_3),
                       parameters::lambda * parameters::lambda * (
                            d_chi_1 * d_theta * sin(parameters::symmetrical::beta[0]) +
                            d_chi_2 * d_theta * sin(parameters::symmetrical::beta[1]) +
                            d_chi_3 * d_theta * sin(parameters::symmetrical::beta[2])
                       ) - cos(parameters::symmetrical::beta[0]) * (parameters::c1 * control[0] - parameters::c2 * d_chi_1) -
                           cos(parameters::symmetrical::beta[1]) * (parameters::c1 * control[1] - parameters::c2 * d_chi_2) -
                           cos(parameters::symmetrical::beta[2]) * (parameters::c1 * control[2] - parameters::c2 * d_chi_3)
                      };


    gsl_matrix_view A = gsl_matrix_view_array(a_data, 3, 3);
    gsl_vector_view b = gsl_vector_view_array(b_data, 3);
    gsl_vector*     gsx = gsl_vector_alloc(3);

    gsl_permutation* p = gsl_permutation_alloc(3);
    int s;
    gsl_linalg_LU_decomp(&A.matrix, p, &s);
    gsl_linalg_LU_solve(&A.matrix, p, &b.vector, gsx);

    N_1.append(gsx->data[0]);
    N_2.append(gsx->data[1]);
    N_3.append(gsx->data[2]);

    gsl_permutation_free(p);
}

const double program_R = 2;

Vector<2> get_program_nu3_coeff(double thetaN_1, double thetaN_2) {
    Vector<2> res;

//    res[0] = 2 * M_PI * parameters::symmetrical::L * (thetaN_2 - 2 * thetaN_1);
//    res[1] = 2 * M_PI * parameters::symmetrical::L * (2 * thetaN_1 - thetaN_2 / 2);
    res[0] = 0;
    res[1] = 2 * M_PI * parameters::symmetrical::L * thetaN_1;

    return res;
}

Vector<3> get_expected_nu(double t, double thetaN_1, double thetaN_2) {
    Vector<3> res;

    auto nu_coeff = get_program_nu3_coeff(thetaN_1, thetaN_2);

    res[0] = 0;
    res[2] = nu_coeff[0] * t + nu_coeff[1];
    res[1] = program_R * res[2] / parameters::symmetrical::L;

    return res;
}

Vector<3> get_expected_dnu(double t, double thetaN_1, double thetaN_2) {
    Vector<3> res;

    auto nu_coeff = get_program_nu3_coeff(thetaN_1, thetaN_2);

    res[0] = 0;
    res[2] = nu_coeff[0];
    res[1] = program_R * res[2] / parameters::symmetrical::L;

    return res;
}

Vector<3> get_program_control_for_symm(double t, double thetaN_1, double thetaN_2) {
    auto nu_v = get_expected_nu(t, thetaN_1, thetaN_2);
    auto dnu_v = get_expected_dnu(t, thetaN_1, thetaN_2);

    double d_1 = parameters::symmetrical::delta[0];
    double d_2 = parameters::symmetrical::delta[1];
    double d_3 = parameters::symmetrical::delta[2];

    double D = parameters::symmetrical::Delta;

    double a_1 = parameters::symmetrical::alpha[0];
    double a_2 = parameters::symmetrical::alpha[1];
    double a_3 = parameters::symmetrical::alpha[2];
    double b_1 = parameters::symmetrical::beta[0];
    double b_2 = parameters::symmetrical::beta[1];
    double b_3 = parameters::symmetrical::beta[2];

    double l = parameters::lambda;
    double L = parameters::symmetrical::L;

    double s_11 = -sin(b_1);
    double s_21 = -sin(b_2);
    double s_31 = -sin(b_3);
    double s_12 = cos(b_1);
    double s_22 = cos(b_2);
    double s_32 = cos(b_3);
    double s_13 = (- D * sin(b_1) + d_1 * cos(a_1 - b_1)) / L;
    double s_23 = (- D * sin(b_2) + d_2 * cos(a_2 - b_2)) / L;
    double s_33 = (- D * sin(b_3) + d_3 * cos(a_3 - b_3)) / L;

    double Xi_raw[9] = {s_11, s_12, s_13,
                        s_21, s_22, s_23,
                        s_31, s_32, s_33};
    gsl_matrix_view Xi = gsl_matrix_view_array(Xi_raw, 3, 3);

    gsl_matrix* dnu_coef = gsl_matrix_alloc(3, 3);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                    l*l, &Xi.matrix, &Xi.matrix,
                    0.0, dnu_coef);

    double eye_3_raw[9] = {1, 0, 0,
                           0, 1, 0,
                           0, 0, 1};
    gsl_matrix_view eye_3 = gsl_matrix_view_array(eye_3_raw, 3, 3);

    gsl_matrix_add(dnu_coef, &eye_3.matrix);


//    qDebug() << dnu_coef->data[0] << " " << dnu_coef->data[1] << " " << dnu_coef->data[2] << '\n';
//    qDebug() << dnu_coef->data[3] << " " << dnu_coef->data[4] << " " << dnu_coef->data[5] << '\n';
//    qDebug() << dnu_coef->data[6] << " " << dnu_coef->data[7] << " " << dnu_coef->data[8] << '\n';

    double dnu_raw[3] = {dnu_v[0], dnu_v[1], dnu_v[2]};
    gsl_vector_view dnu = gsl_vector_view_array(dnu_raw, 3);

    gsl_vector* rp_dnu = gsl_vector_alloc(3);
    gsl_blas_dgemv(CblasNoTrans, 1.0, dnu_coef, &dnu.vector, 0, rp_dnu);

    gsl_matrix* nu_coef = gsl_matrix_alloc(3, 3);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                    parameters::c2, &Xi.matrix, &Xi.matrix,
                    0.0, nu_coef);

    double nu_raw[3] = {nu_v[0], nu_v[1], nu_v[2]};
    gsl_vector_view nu = gsl_vector_view_array(nu_raw, 3);

    gsl_vector* rp_nu = gsl_vector_alloc(3);
    gsl_blas_dgemv(CblasNoTrans, 1.0, nu_coef, &nu.vector, 0, rp_nu);

    double rp_numixed_raw[3] = {- nu_v[1] * nu_v[2] / L, nu_v[0] * nu_v[2] / L, 0};
    gsl_vector_view rp_numixed = gsl_vector_view_array(rp_numixed_raw, 3);

    gsl_vector* rp = rp_nu;
    gsl_vector_add(rp, rp_dnu);
    gsl_vector_add(rp, &rp_numixed.vector);

    gsl_matrix* lp = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                    parameters::c1, &Xi.matrix, &eye_3.matrix,
                    0.0, lp);

    gsl_vector *x = gsl_vector_alloc(3);

    int s;

    gsl_permutation *p = gsl_permutation_alloc(3);

    gsl_linalg_LU_decomp(lp, p, &s);

    gsl_linalg_LU_solve(lp, p, rp, x);

    Vector<3> res = {x->data[0], x->data[1], x->data[2]};

    qDebug() << res[0] << " " << res[1] << " " << res[2] << '\n';

    gsl_permutation_free(p);

    return res;
}

Vector<6> DOPRI8_symmetrical_plot(double t_left, double t_right,
                    double thetaN_1, double thetaN_2,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec,
                    QVector<double> &P_real, QVector<double> &P_advice,
                    QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3)
{
    double h = (t_right - t_left) / 1e7;
    double h_new;
    bool last_flag = false;

    double tl = t_left;
    Vector<6> xl;

    Vector<6> stepx;

    Vector<6> errx;
    double err, coef;
    double coefmax = 5;

    Vector<6> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

    auto initial_nu = get_expected_nu(tl, thetaN_1, thetaN_2);
    Vector<3> control = get_program_control_for_symm(tl, thetaN_1, thetaN_2);
    xl = {initial_nu[0], initial_nu[1], initial_nu[2], 0, 0, 0};

    t_vec.append(tl);
    nu1_vec.append(xl[0]);
    nu2_vec.append(xl[1]);
    nu3_vec.append(xl[2]);
    x_vec.append(xl[3]);
    y_vec.append(xl[4]);
    theta_vec.append(xl[5]);

    compute_P(tl, xl, control, P_real, P_advice);
    compute_N(tl, xl, control, N_1, N_2, N_3);

    while (tl + h < t_right || last_flag)
    {
        qDebug() << "nu1: " << xl[0] << "nu2: " << xl[1] << "nu3: " << xl[2] << "\n";

        control = get_program_control_for_symm(tl, thetaN_1, thetaN_2);

        k1 = rightpart(tl, xl, control);

        k2 = rightpart(tl + h / 18, xl + h * k1 / 18, control);

        k3 = rightpart(tl + h / 12, xl + h * (k1  / 48 + k2 / 16) , control);

        k4 = rightpart(tl + h / 8, xl + h * (k1  / 32 + k3 * 3 / 32) , control);

        k5 = rightpart(tl + h * 5 / 16, xl + h * (k1  * 5 / 16 + k3 * (-75) / 64 + k4 * 75 / 64) , control);

        k6 = rightpart(tl + h * 3 / 8, xl + h * (k1  * 3 / 80 + k4 * 3 / 16 + k5 * 3 / 20) , control);

        k7 = rightpart(tl + h * 59. / 400., xl + h * (k1  * (29443841. / 614563906.) + k4 * (77736538. / 692538347.) + k5 * ((-28693883.) / 1125000000.) + 
                k6 * (23124283. / 1800000000.) ), control);

        k8 = rightpart(tl + h * 93. / 200., xl + h * (k1  * (16016141. / 946692911.) + k4 * (61564180. / 158732637.) + k5 * (22789713. / 633445777.) + 
                k6 * (545815736. / 2771057229.) + k7 * ((-180193667.) / 1043307555.) ), control);

        k9 = rightpart(tl + h * (5490023248. / 9719169821.), xl + h * (k1  * (39632708. / 573591083.) + k4 * ((-433636366.) / 683701615.) + k5 * ((-421739975.) / 2616292301.) +
                k6 * (100302831. / 723423059.) + k7 * (790204164. / 839813087.) + k8 * (800635310. / 3783071287.) ), control);

        k10 = rightpart(tl + h * 13 / 20, xl + h * (k1  * (246121993. / 1340847787.) + k4 * ((-37695042795.) / 15268766246.) + k5 * ((-309121744.) / 1061227803.) +
                k6 * ((-12992083.) / 490766935.) + k7 * (6005943493. / 2108947869.) + k8 * (393006217. / 1396673457) + k9 * (123872331. / 1001029789.) ), control);

        k11 = rightpart(tl + h * (1201146811. / 1299019798.), xl + h * (k1  * ((-1028468189.) / 846180014.) + k4 * (8478235783. / 508512852.) + k5 * (1311729495. / 1432422823.) +
                k6 * ((-10304129995.) / 1701304382.) + k7 * ((-48777925059.) / 3047939560.) + k8 * (15336726248. / 1032824649.) + k9 * ((-45442868181.) / 3398467696.) + 
                k10 * (3065993473. / 597172653.) ), control);

        k12 = rightpart(tl + h, xl + h * (k1  * (185892177. / 718116043.) + k4 * ((-3185094517.) / 667107341.) + k5 * ((-477755414.) / 1098053517.) + 
                 k6 * ((-703635378.) / 230739211.) + k7 * (5731566787. / 1027545527.) + k8 * (5232866602. / 850066563.) + k9 * ((-4093664535.) / 808688257.) + 
                k10 * (3962137247. / 1805957418.) + k11 * (65686358. / 487910083.) ), control);

        k13 = rightpart(tl + h, xl + h * (k1  * (403863854. / 491063109.) + k4 * ((-5068492393.) / 434740067.) + k5 * ((-411421997.) / 543043805.) + 
                 k6 * (652783627. / 914296604.) + k7 * (11173962825. / 925320556.) + k8 * ((-13158990841.) / 6184727034.) + k9 * (3936647629. / 1978049680.) + 
                k10 * ((-160528059.) / 685178525.) + k11 * (248638103. / 1413531060.) ), control);

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
            nu1_vec.append(xl[0]);
            nu2_vec.append(xl[1]);
            nu3_vec.append(xl[2]);
            x_vec.append(xl[3]);
            y_vec.append(xl[4]);
            theta_vec.append(xl[5]);
            compute_P(tl, xl, control, P_real, P_advice);
            compute_N(tl, xl, control, N_1, N_2, N_3);

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
                nu1_vec.append(xl[0]);
                nu2_vec.append(xl[1]);
                nu3_vec.append(xl[2]);
                x_vec.append(xl[3]);
                y_vec.append(xl[4]);
                theta_vec.append(xl[5]);
                compute_P(tl, xl, control, P_real, P_advice);
                compute_N(tl, xl, control, N_1, N_2, N_3);

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

            qDebug() << "last" << '\n';
        }
        else
            last_flag = false;
    }

    return xl;
}
