#include "parameters.h"
#include "vector.h"
#include <QVector>
#include <QDebug>
#include <math.h>
#include <gsl/gsl_linalg.h>

//nu1, nu2, nu3, x, y, theta

//Vector<3> compute_W(double t) {
//    double eps_0 = 0.01;
//    double omega_0 = 0.1;

//    double R = 2;

//    double L = parameters::symmetrical::Lambda;
//    double r = parameters::symmetrical::rho;
//    double c2 = parameters::c2;
//    double c1 = parameters::c1;
//    double l = parameters::lambda;

//    Vector<3> W;

//    qDebug() << M_PI;

//    W[0] = - 8 * M_PI * M_PI * R / (3 * l * l + 2) * (eps_0 * t + omega_0) * (eps_0 * t + omega_0);
//    W[1] = 2 * M_PI * R * eps_0 + 3 * c2 / (3 * l * l + 2) * 2 * M_PI * R * (eps_0 * t + omega_0);
//    W[2] = 2 * M_PI * L * eps_0 + 3 * c2 * r * r / (3 * l * l * r * r + L * L) * 2 * M_PI * L * (eps_0 * t + omega_0);

//    return W;
//}

//Vector<3> compute_U(double t) {
//    double eps_0 = 0.01;
//    double omega_0 = 0.1;

//    double R = 2;

//    double L = parameters::symmetrical::Lambda;
//    double r = parameters::symmetrical::rho;
//    double c2 = parameters::c2;
//    double c1 = parameters::c1;
//    double l = parameters::lambda;

//    Vector<3> W;
//    W = compute_W(t);

//    double w0 = W[0] / c1 * (3 * l * l + 2);
//    double w1 = W[1] / c1 / sqrt(3) * (3 * l * l + 2);
//    double w2 = W[2] / c1 / r / L * (3 * l * l * r * r + L * L);

//    Vector<3> U;
//    U[0] = w0 / 6 + w1 / 2 + w2 / 3;
//    U[1] = - w0 / 3 + w2 / 3;
//    U[2] = w0 / 6 - w1 / 2 + w2 / 3;

//    return U;
//}

Vector<3> get_desired_control(double dnu_1, double dnu_2, double dnu_3,
                           double nu_1, double nu_2, double nu_3)
{
    double L = parameters2::Lambda;
    double l = parameters2::lambda;

    double m = parameters2::m;
    double h = parameters2::h;
    double R = parameters2::R;
    double rho = parameters2::rho;

    double c_1 = parameters2::c_1;
    double c_2 = parameters2::c_2;

    double w1 = (2 * m * R * R + 3 * l * l) / (c_1 * R) * (dnu_1 - 2 * m * R * R / (2 * m * R * R + 3 * l * l) * nu_2 * nu_3 / L + \
                                                           3 * c_2 / (2 * m * R * R + 3 * l * l) * nu_1);
    double w2 = (2 * m * R * R + 3 * l * l) / (sqrt(3.0) * c_1 * R) * (dnu_2 + 2 * m * R * R / (2 * m * R * R + 3 * l * l) * nu_1 * nu_3 / L + \
                                                                       3 * c_2 / (2 * m * R * R + 3 * l * l) * nu_2);
    double w3 = (L * L * R * R + 3 * l * l * rho * rho) / (c_1 * L * rho * R) * (dnu_3 + 3 * c_2 * rho * rho / (L * L * R * R + 3 * l * l * rho * rho) * nu_3);

    double U1 = w1 / 6 + w2 / 2 + w3 / 3;
    double U2 = - w1 / 3 + w3 / 3;
    double U3 = w1 / 6 - w2 / 2 + w3 / 3;

    Vector<3> U;
        U[0] = U1;
        U[1] = U2;
        U[2] = U3;

        return U;
}

Vector<6> get_dangerous(double t) {
    double dH = 50;
    double H = dH * (t + 1e-13);

    double L = parameters2::Lambda;
    double l = parameters2::lambda;

    double m = parameters2::m;
    double h = parameters2::h;
    double R = parameters2::R;
    double rho = parameters2::rho;


    double c_1 = parameters2::c_1;
    double c_2 = parameters2::c_2;

    double dnu_3 = 0.03;
    double nu_3 = dnu_3 * t;

    double k = sqrt(2.0) * l * l * L / 2.0 * sqrt(L * L * R * R + 3 * l * l * rho * rho) / sqrt(2 * m * R * R + 3 * l * l);

    double nu_1 = - l * l * nu_3 / 2.0 / k;
    double dnu_1 = - l * l * dnu_3 / 2.0 / k;

    double nu_2 = - sqrt(3.0) * l * l * nu_3 / 2.0 / k;
    double dnu_2 = - sqrt(3.0) * l * l * dnu_3 / 2.0 / k;

    Vector<6> res;
    res[0] = dnu_1;
    res[1] = dnu_2;
    res[2] = dnu_3;
    res[3] = nu_1;
    res[4] = nu_2;
    res[5] = nu_3;

    qDebug() << res[0] << '\n';
    return res;
}

Vector<3> compute_U(double t) {
    double eps_0 = 0.04;
    double R_traj = 1;

    double L = parameters2::Lambda;
    double l = parameters2::lambda;

    double m = parameters2::m;
    double h = parameters2::h;
    double R = parameters2::R;
    double rho = parameters2::rho;


    double c_1 = parameters2::c_1;
    double c_2 = parameters2::c_2;

    double dnu_3 = 2 * M_PI * L * eps_0;
    double nu_3 = dnu_3 * t;
    double dnu_2 = R_traj / L * dnu_3;
    double nu_2 = R_traj / L * nu_3;

//    qDebug() << "U1 " << U1 << "U2 " << U2 << "U3 " << U3 << '\n';

    //    Vector<3> W;
    //    W = compute_W(t);

    //    double w0 = W[0] / c1 * (3 * l * l + 2);
    //    double w1 = W[1] / c1 / sqrt(3) * (3 * l * l + 2);
    //    double w2 = W[2] / c1 / r / L * (3 * l * l * r * r + L * L);

    //    Vector<3> U;
    //    U[0] = w0 / 6 + w1 / 2 + w2 / 3;
    //    U[1] = - w0 / 3 + w2 / 3;
    //    U[2] = w0 / 6 - w1 / 2 + w2 / 3;

    //    return U;

//    Vector<3> U;
//    U[0] = U1;
//    U[1] = U2;
//    U[2] = U3;
//    U[0] = -48;
//    U[1] = 0;
//    U[2] = 48;

//    auto dang = get_dangerous(t);
    return get_desired_control(0, dnu_2, dnu_3, 0, nu_2, nu_3);
//    return {-20, 30, -12};
}


Vector<6> rightpart(double t, Vector<6> x)
{
    Vector<6> res;

    Vector<3> control;
    control = compute_U(t);

    double U_1 = control[0];
    double U_2 = control[1];
    double U_3 = control[2];

    double L = parameters2::Lambda;
    double l = parameters2::lambda;

    double m = parameters2::m;
    double h = parameters2::h;
    double R = parameters2::R;
    double rho = parameters2::rho;


    double c_1 = parameters2::c_1;
    double c_2 = parameters2::c_2;

//    rez[0] = (x[1] * x[2] / parameters::symmetrical::Lambda + parameters::c1 * (u[0] / 2 - u[1] + u[2] / 2)
//            - 3 * parameters::c2 * x[0] / 2) / parameters::symmetrical::A1;
//    rez[0] = 0;
//    rez[1] = (- x[0] * x[2] / parameters::symmetrical::Lambda + parameters::c1 * sqrt(3) * (u[0]  - u[2]) / 2
//            - 3 * parameters::c2 * x[1] / 2) / parameters::symmetrical::A2;
//    rez[2] = (parameters::c1 * parameters::symmetrical::rho * (u[0] + u[1] + u[2]) / parameters::symmetrical::Lambda
//            - 3 * parameters::c2 * parameters::symmetrical::rho * parameters::symmetrical::rho * x[2]
//            / (parameters::symmetrical::Lambda * parameters::symmetrical::Lambda)) / parameters::symmetrical::A3;
    res[0] = 2 * m * R * R / (2 * m * R * R + 3 * l * l) * x[1] * x[2] / L - 3 * c_2 / (2 * m * R * R + 3 * l * l) * x[0] + c_1 * R / (2 * m * R * R + 3 * l * l) * (U_1 - 2 * U_2 + U_3);
//    rez[0] = 0;
    res[1] = - 2 * m * R * R / (2 * m * R * R + 3 * l * l) * x[0] * x[2] / L - 3 * c_2 / (2 * m * R * R + 3 * l * l) * x[1] + sqrt(3.0) * c_1 * R / (2 * m * R * R + 3 * l * l) * (U_1 - U_3);
    res[2] = - 3 * c_2 * rho * rho / (L * L * R * R + 3 * l * l * rho * rho) * x[2] + c_1 * L * rho * R / (L * L * R * R + 3 * l * l * rho * rho) * (U_1 + U_2 + U_3);
    res[3] = x[0] * cos(x[5]) - x[1] * sin(x[5]);
    res[4] = x[0] * sin(x[5]) + x[1] * cos(x[5]);
    res[5] = x[2] / L;

    return res;
}

void compute_P(double t, Vector<6> x, QVector<double> &P_real, QVector<double> &P_advice,
               QVector<double> &PT_advice)
{
    Vector<3> control;
    control = compute_U(t);

    double U_1 = control[0];
    double U_2 = control[1];
    double U_3 = control[2];

    double L = parameters2::Lambda;
    double l = parameters2::lambda;

    double m = parameters2::m;
    double h = parameters2::h;
    double R = parameters2::R;
    double rho = parameters2::rho;


    double c_1 = parameters2::c_1;
    double c_2 = parameters2::c_2;

    double g = 9.81;

    double fru_0 = sqrt(3.0) * m * rho * g / 3.0 / c_1 * (2 * m * R * R + 3 * l * l) / (2 * m * h * R + 3 * l * l);
    double fru_1 = sqrt(3.0) * c_2 / c_1 / R;
    double fru_2 = sqrt(3.0) * l * l / c_1 / L / R;


//    double H = (1 + 3.0 * l * l / 2) * (x[0] * x[0] + x[1] * x[1]) + (1 + 3.0 * r * r * l * l / L / L) * x[2] * x[2];


//    double good_advice = fru_0 / sqrt(2.0) - sqrt(3.0) * r / c1 * \
//                         sqrt(3 * l * l / r / r + L * L) / sqrt(3 * l * l + 2) * \
//                         (l * l / (3 * l * l / r / r + L * L) * H + c2 * c2 / 2 / l / l);

//    double l = parameters::lambda;
//    double good_advice = fru_0 / sqrt(2.0) - sqrt(2.0) * fru_1 * sqrt(H / (3.0 * l * l + 2.0));
//    qDebug() << "H" << H << '\n';
//    qDebug() << "good_advice" << good_advice << '\n';

    P_real.append(sqrt(U_1 * U_1 + U_2 * U_2 + U_3 * U_3));

    double v_s = sqrt(x[0] * x[0] + x[1] * x[1]);
    double v_advice = fru_0 / sqrt(2.0) - v_s / sqrt(2.0) * sqrt(fru_2 * fru_2 * x[2] * x[2] + fru_1 * fru_1);
    P_advice.append(std::max(v_advice, 0.0));

    double H = (m + 3.0 * l * l / 2 / R / R) * (x[0] * x[0] + x[1] * x[1]) + (1 + 3.0 * rho * rho * l * l / L / L / R / R) * x[2] * x[2];
    double T_advice = fru_0 / sqrt(2.0) - sqrt(3.0) / c_1 / R * \
                             sqrt(L * L * R * R + 3 * l * l * rho * rho) / sqrt(2 * m * R * R + 3 * l * l) * \
                            (l * l * R * R/ (3 * l * l * rho * rho + L * L * R * R) * H + c_2 * c_2 / 2 / l / l);
    PT_advice.append(std::max(T_advice, 0.0));
}

void compute_N(double t, Vector<6> x,
               QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3)
{
    Vector<3> control;

    control = compute_U(t);

    double U_1 = control[0];
    double U_2 = control[1];
    double U_3 = control[2];

    double L = parameters2::Lambda;
    double l = parameters2::lambda;

    double m = parameters2::m;
    double h = parameters2::h;
    double R = parameters2::R;
    double rho = parameters2::rho;


    double c_1 = parameters2::c_1;
    double c_2 = parameters2::c_2;

    double nu_1 = x[0];
    double nu_2 = x[1];
    double nu_3 = x[2];

    double g = parameters2::g;

    N_1.append(
                m * g / 3 + (2 * m * h * R + 3 * l * l) / rho / (2 * m * R * R + 3 * l * l) * \
                            ( \
                                l * l / L / R * (nu_1 + sqrt(3.0) * nu_2) / 2 * nu_3 + \
                                c_2 / R * (sqrt(3.0) * nu_1 - nu_2) / 2 + \
                                c_1 * (U_2 - U_3) / sqrt(3.0) \
                            ) \
               );
    N_2.append(
                m * g / 3 + (2 * m * h * R + 3 * l * l) / rho / (2 * m * R * R + 3 * l * l) * \
                            ( \
                                - l * l / L / R * nu_1 * nu_3 + \
                                c_2 / R * nu_2 + \
                                c_1 * (U_3 - U_1) / sqrt(3.0) \
                            ) \
               );
    N_3.append(
                m * g / 3 + (2 * m * h * R + 3 * l * l) / rho / (2 * m * R * R + 3 * l * l) * \
                            ( \
                                l * l / L / R * (nu_1 - sqrt(3.0) * nu_2) / 2 * nu_3 + \
                                c_2 / R * (- sqrt(3.0) * nu_1 - nu_2) / 2 + \
                                c_1 * (U_1 - U_2) / sqrt(3.0) \
                            ) \
               );
}

Vector<6> DOPRI8_symmetrical_plot(double t_left, double t_right, Vector<6> initial_values,
                    Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec,
                    QVector<double> &P_real, QVector<double> &P_advice, QVector<double> &PT_advice,
                    QVector<double> &N_1, QVector<double> &N_2, QVector<double> &N_3,
                    QVector<double> &U1, QVector<double> &U2, QVector<double> &U3)
{
    double h = (t_right - t_left) / 1e7;
    double h_new;
    bool switch_flag = false;
    bool last_flag = false;

    double tl = t_left;
    Vector<6> xl = initial_values;


//    double L = parameters::symmetrical::Lambda;
//    double r = parameters::symmetrical::rho;
//    double c2 = parameters::c2;
//    double c1 = parameters::c1;
//    double l = parameters::lambda;


//    double eps_0 = 0.01;
//    double omega_0 = 0.1;
//    double R = 2;

//    xl[0] = 0;
//    xl[2] = 2 * M_PI * L * omega_0;
//    xl[1] = R / L * xl[2];
    xl[0] = 0;
    xl[1] = 0;
    xl[2] = 0;

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
    compute_P(tl, xl, P_real, P_advice, PT_advice);
    compute_N(tl, xl, N_1, N_2, N_3);

    Vector<3> u;
    u = compute_U(tl);
    U1.append(u[0]);
    U2.append(u[1]);
    U3.append(u[2]);

    while (tl + h < t_right || last_flag)
    {

        //switch point
        if (tl < t_sw && tl + h > t_sw && !switch_flag)
        {
            h = t_sw - tl;
            switch_flag = true;
        }

        k1 = rightpart(tl, xl);

        k2 = rightpart(tl + h / 18, xl + h * k1 / 18);

        k3 = rightpart(tl + h / 12, xl + h * (k1  / 48 + k2 / 16));

        k4 = rightpart(tl + h / 8, xl + h * (k1  / 32 + k3 * 3 / 32));

        k5 = rightpart(tl + h * 5 / 16, xl + h * (k1  * 5 / 16 + k3 * (-75) / 64 + k4 * 75 / 64));

        k6 = rightpart(tl + h * 3 / 8, xl + h * (k1  * 3 / 80 + k4 * 3 / 16 + k5 * 3 / 20));

        k7 = rightpart(tl + h * 59. / 400., xl + h * (k1  * (29443841. / 614563906.) + k4 * (77736538. / 692538347.) + k5 * ((-28693883.) / 1125000000.) + 
                k6 * (23124283. / 1800000000.) ));

        k8 = rightpart(tl + h * 93. / 200., xl + h * (k1  * (16016141. / 946692911.) + k4 * (61564180. / 158732637.) + k5 * (22789713. / 633445777.) + 
                k6 * (545815736. / 2771057229.) + k7 * ((-180193667.) / 1043307555.) ));

        k9 = rightpart(tl + h * (5490023248. / 9719169821.), xl + h * (k1  * (39632708. / 573591083.) + k4 * ((-433636366.) / 683701615.) + k5 * ((-421739975.) / 2616292301.) +
                k6 * (100302831. / 723423059.) + k7 * (790204164. / 839813087.) + k8 * (800635310. / 3783071287.) ));

        k10 = rightpart(tl + h * 13 / 20, xl + h * (k1  * (246121993. / 1340847787.) + k4 * ((-37695042795.) / 15268766246.) + k5 * ((-309121744.) / 1061227803.) +
                k6 * ((-12992083.) / 490766935.) + k7 * (6005943493. / 2108947869.) + k8 * (393006217. / 1396673457) + k9 * (123872331. / 1001029789.) ));

        k11 = rightpart(tl + h * (1201146811. / 1299019798.), xl + h * (k1  * ((-1028468189.) / 846180014.) + k4 * (8478235783. / 508512852.) + k5 * (1311729495. / 1432422823.) +
                k6 * ((-10304129995.) / 1701304382.) + k7 * ((-48777925059.) / 3047939560.) + k8 * (15336726248. / 1032824649.) + k9 * ((-45442868181.) / 3398467696.) + 
                k10 * (3065993473. / 597172653.) ));

        k12 = rightpart(tl + h, xl + h * (k1  * (185892177. / 718116043.) + k4 * ((-3185094517.) / 667107341.) + k5 * ((-477755414.) / 1098053517.) + 
                 k6 * ((-703635378.) / 230739211.) + k7 * (5731566787. / 1027545527.) + k8 * (5232866602. / 850066563.) + k9 * ((-4093664535.) / 808688257.) + 
                k10 * (3962137247. / 1805957418.) + k11 * (65686358. / 487910083.) ));

        k13 = rightpart(tl + h, xl + h * (k1  * (403863854. / 491063109.) + k4 * ((-5068492393.) / 434740067.) + k5 * ((-411421997.) / 543043805.) + 
                 k6 * (652783627. / 914296604.) + k7 * (11173962825. / 925320556.) + k8 * ((-13158990841.) / 6184727034.) + k9 * (3936647629. / 1978049680.) + 
                k10 * ((-160528059.) / 685178525.) + k11 * (248638103. / 1413531060.) ));

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
            compute_P(tl, xl, P_real, P_advice, PT_advice);
            compute_N(tl, xl, N_1, N_2, N_3);
            u = compute_U(tl);
            U1.append(u[0]);
            U2.append(u[1]);
            U3.append(u[2]);

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
                compute_P(tl, xl, P_real, P_advice, PT_advice);
                compute_N(tl, xl, N_1, N_2, N_3);
                u = compute_U(tl);
                U1.append(u[0]);
                U2.append(u[1]);
                U3.append(u[2]);

                coefmax = 5;
            }
            else
                coefmax = 1;

            h = h_new;
        }

        //for correct plotting; there could be too few points
        if (h > (t_right - t_left) / 500)
        {
            h = (t_right - t_left) / 500;
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
