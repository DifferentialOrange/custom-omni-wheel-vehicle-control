#include "parameters.h"
#include "vector.h"
#include <QVector>
#include <QDebug>
#include <math.h>
#include <gsl/gsl_linalg.h>

//nu1, nu2, nu3, x, y, theta

Vector<6> rightpart(double t, Vector<6> x, Vector<3> control_minus, Vector<3> control_plus, double t_sw)
{
    Vector<6> rez;
    Vector<3> u;


    if (t < t_sw)
        u = control_minus;
    else
        u = control_plus;


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

void compute_P(double t, Vector<6> x, Vector<3> control_minus, Vector<3> control_plus, double t_sw, QVector<double> &P_real, QVector<double> &P_advice)
{
    Vector<3> u;

    if (t < t_sw)
        u = control_minus;
    else
        u = control_plus;

    P_real.append(sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));

    double g = 9.81;

//    double h = (1 + 3.0 * parameters::lambda * parameters::lambda / 2) * (x[0] * x[0] + x[1] * x[1]) +
//            (1 + 3.0 * parameters::symmetrical::rho * parameters::symmetrical::rho * parameters::lambda * parameters::lambda /
//             parameters::symmetrical::Lambda / parameters::symmetrical::Lambda) * x[2] * x[2];

    double v_S = sqrt(x[0] * x[0] + x[1] * x[1]);

//    double good_advice = g * parameters::symmetrical::rho / parameters::c1 / sqrt(6) -
//                (3 + sqrt(3)) * (parameters::lambda * parameters::lambda * h / sqrt(parameters::symmetrical::Lambda * parameters::symmetrical::Lambda +
//                                                                                    3 * parameters::symmetrical::rho * parameters::symmetrical::rho *
//                                                                                    parameters::lambda * parameters::lambda) +
//                                 parameters::c2 * sqrt(h)) / 2 / parameters::c1 / sqrt(2 + 3 * parameters::lambda * parameters::lambda);
    double good_advice = sqrt(6) * g * parameters::symmetrical::rho * (3 * parameters::lambda * parameters::lambda + 2) /
                            (6 * parameters::c1 * (3 * parameters::lambda * parameters::lambda + 2 * 1)) -
                            sqrt(6) * sqrt(parameters::lambda * parameters::lambda * parameters::lambda * parameters::lambda *
                                           x[2] * x[2] + parameters::c2 * parameters::c2) * v_S / 2 / parameters::c1;


    P_advice.append(good_advice);
}

void compute_N(double t, Vector<6> x, Vector<3> control_minus, Vector<3> control_plus, double t_sw,
               QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3)
{
    Vector<3> control;

    if (t < t_sw)
        control = control_minus;
    else
        control = control_plus;

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
}

Vector<6> DOPRI8_symmetrical_plot(double t_left, double t_right, Vector<6> initial_values,
                    Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec,
                    QVector<double> &P_real, QVector<double> &P_advice,
                    QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3)
{
    double h = (t_right - t_left) / 1e7;
    double h_new;
    bool switch_flag = false;
    bool last_flag = false;

    double tl = t_left;
    Vector<6> xl = initial_values;

    Vector<6> stepx;

    Vector<6> errx;
    double err, coef;
    double coefmax = 5;

    Vector<6> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

    t_vec.append(tl);
    nu1_vec.append(xl[0]);
    nu2_vec.append(xl[1]);
    nu3_vec.append(xl[2]);
    x_vec.append(xl[3]);
    y_vec.append(xl[4]);
    theta_vec.append(xl[5]);
    compute_P(tl, xl, control_minus, control_plus, t_sw, P_real, P_advice);
    compute_N(tl, xl, control_minus, control_plus, t_sw, N_1, N_2, N_3);

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
            nu1_vec.append(xl[0]);
            nu2_vec.append(xl[1]);
            nu3_vec.append(xl[2]);
            x_vec.append(xl[3]);
            y_vec.append(xl[4]);
            theta_vec.append(xl[5]);
            compute_P(tl, xl, control_minus, control_plus, t_sw, P_real, P_advice);
            compute_N(tl, xl, control_minus, control_plus, t_sw, N_1, N_2, N_3);

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
                compute_P(tl, xl, control_minus, control_plus, t_sw, P_real, P_advice);
                compute_N(tl, xl, control_minus, control_plus, t_sw, N_1, N_2, N_3);

                coefmax = 5;
            }
            else
                coefmax = 1;

            h = h_new;
        }

        //for correct plotting; there could be too few points
        if (h > (t_right - t_left) / 10)
        {
            h = (t_right - t_left) / 10;
        }

        //last step correction
        if (tl + h >= t_right && last_flag == false)
        {
            h = t_right - tl;
            last_flag = true;
        }
        else
            last_flag = false;
    }

    return xl;
}
