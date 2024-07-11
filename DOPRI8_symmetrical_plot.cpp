#include "parameters.h"
#include "vector.h"
#include <QVector>
#include <QDebug>
#include <math.h>
#include <gsl/gsl_linalg.h>

double rho(double chi_2) {
    double r = 0.05;
    double n = 4.0;

    return - 1.0 / r / (cos(chi_2) - cos(M_PI / n));
}

double drho(double chi_2) {
    double r = 0.05;
    double n = 4.0;

    return - sin(chi_2) / r / (cos(chi_2) - cos(M_PI / n)) / (cos(chi_2) - cos(M_PI / n));
}

// nu_1, nu_2, nu_3, x, y, theta, chi_2
Vector<7> rightpart(double t, Vector<7> x, double eps)
{
    Vector<7> res;

    double d_ll = 0.1;
    double d_2 = 0.15;
    double D = 0.05;
//    double h = 0.05;
    double r = 0.05;
    double m = 3;
    double L = sqrt(5e-2);
    double l = sqrt(5e-4);

    double a_1_triag = m + l * l / r / r;
    double a_2_triag = m + 2 * l * l / r / r;
    double a_3_triag = 1 + l * l * (2 * d_ll * d_ll + (d_2 - D) * (d_2 - D)) / L / L / r / r;
    double k = l * l / a_3_triag * (d_2 - D) / L / r / r;

    double nu_1 = x[0];
    double nu_2 = x[1];
    double nu_3 = x[2];
    double theta = x[5];
    double chi_2 = x[6];

    double a_11 = a_1_triag;
    double a_12 = 0;
    double a_13 = - k * a_3_triag;
    double a_21 = a_12;
    double a_22 = a_2_triag + eps * eps * rho(chi_2) * rho(chi_2);
    double a_23 = eps * eps * rho(chi_2) * sin(chi_2) / L;
    double a_31 = a_13;
    double a_32 = a_23;
    double a_33 = a_3_triag;

    double b_1 = (m / L - eps * eps * rho(chi_2) * cos(chi_2) / L / r) * nu_2 * nu_3;
    double b_2 = - (m / L - eps * eps * rho(chi_2) * cos(chi_2) / L / r) * nu_1 * nu_3 + \
                 eps * eps * rho(chi_2) * drho(chi_2) / r * nu_1 * nu_2 - \
                 eps * eps * (d_2 - D) * rho(chi_2) * drho(chi_2) / L / r * nu_2 * nu_3 - \
                 eps * eps * (d_2 - D) * rho(chi_2) * cos(chi_2) / L / L / r * nu_3 * nu_3;
    double b_3 = eps * eps * (rho(chi_2) * cos(chi_2) + drho(chi_2) * sin(chi_2)) / L / r * nu_1 * nu_2 - \
                 eps * eps * (d_2 - D) * drho(chi_2) * sin(chi_2) / L / L / r * nu_2 * nu_3;

    int s;
    // a y = b
    gsl_matrix * a = gsl_matrix_alloc (3, 3);
    gsl_matrix_set (a, 0, 0, a_11);
    gsl_matrix_set (a, 0, 1, a_12);
    gsl_matrix_set (a, 0, 2, a_13);
    gsl_matrix_set (a, 1, 0, a_21);
    gsl_matrix_set (a, 1, 1, a_22);
    gsl_matrix_set (a, 1, 2, a_23);
    gsl_matrix_set (a, 2, 0, a_31);
    gsl_matrix_set (a, 2, 1, a_32);
    gsl_matrix_set (a, 2, 2, a_33);
    gsl_vector * y = gsl_vector_alloc (3);
    gsl_vector * b = gsl_vector_alloc (3);
    gsl_vector_set (b, 0, b_1);
    gsl_vector_set (b, 1, b_2);
    gsl_vector_set (b, 2, b_3);
    gsl_permutation * p = gsl_permutation_alloc (3);
    gsl_linalg_LU_decomp (a, p, &s);
    gsl_linalg_LU_solve (a, p, b, y);

    for (int i = 0; i < 3; i++) {
        res[i] = gsl_vector_get(y, i);
    }
    res[3] = nu_1 * cos(theta) - nu_2 * sin(theta);
    res[4] = nu_1 * sin(theta) + nu_2 * cos(theta);
    res[5] = nu_3 / L;
    res[6] = - 1.0 / r * nu_1 + (d_2 - D) / L / r * nu_3;

    gsl_permutation_free (p);
    gsl_matrix_free (a);
    gsl_vector_free (y);
    gsl_vector_free (b);

    return res;
}

void DOPRI8_symmetrical_plot(double t_left, double t_right,
                    QVector<double> &t_vec, QVector<double> &nu1_vec, QVector<double> &nu2_vec,
                    QVector<double> &nu3_vec, QVector<double> &x_vec, QVector<double> &y_vec,
                    QVector<double> &theta_vec, QVector<double> &chi_2_vec, double eps)
{
    double h = (t_right - t_left) / 1e7;
    double h_new;
    bool last_flag = false;

    double tl = t_left;
    Vector<7> xl;

    // nu_1, nu_2, nu_3, x, y, theta, chi_2
    xl[0] = 0.005;
    xl[1] = - 0.3;
    xl[2] = 0.005 * 3.0 / 2.0 * sqrt(5e-2);
    xl[3] = 0;
    xl[4] = 0;
    xl[5] = 0;
    xl[6] = 0;

    Vector<7> stepx;

    Vector<7> errx;
    double err, coef;
    double coefmax = 5;

    Vector<7> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

    t_vec.append(tl);
    nu1_vec.append(xl[0]);
    nu2_vec.append(xl[1]);
    nu3_vec.append(xl[2]);
    x_vec.append(xl[3]);
    y_vec.append(xl[4]);
    theta_vec.append(xl[5]);
    chi_2_vec.append(xl[6]);

    while (tl + h < t_right || last_flag)
    {
        k1 = rightpart(tl, xl, eps);

        k2 = rightpart(tl + h / 18, xl + h * k1 / 18, eps);

        k3 = rightpart(tl + h / 12, xl + h * (k1  / 48 + k2 / 16), eps);

        k4 = rightpart(tl + h / 8, xl + h * (k1  / 32 + k3 * 3 / 32), eps);

        k5 = rightpart(tl + h * 5 / 16, xl + h * (k1  * 5 / 16 + k3 * (-75) / 64 + k4 * 75 / 64), eps);

        k6 = rightpart(tl + h * 3 / 8, xl + h * (k1  * 3 / 80 + k4 * 3 / 16 + k5 * 3 / 20), eps);

        k7 = rightpart(tl + h * 59. / 400., xl + h * (k1  * (29443841. / 614563906.) + k4 * (77736538. / 692538347.) + k5 * ((-28693883.) / 1125000000.) + 
                k6 * (23124283. / 1800000000.) ), eps);

        k8 = rightpart(tl + h * 93. / 200., xl + h * (k1  * (16016141. / 946692911.) + k4 * (61564180. / 158732637.) + k5 * (22789713. / 633445777.) + 
                k6 * (545815736. / 2771057229.) + k7 * ((-180193667.) / 1043307555.) ), eps);

        k9 = rightpart(tl + h * (5490023248. / 9719169821.), xl + h * (k1  * (39632708. / 573591083.) + k4 * ((-433636366.) / 683701615.) + k5 * ((-421739975.) / 2616292301.) +
                k6 * (100302831. / 723423059.) + k7 * (790204164. / 839813087.) + k8 * (800635310. / 3783071287.) ), eps);

        k10 = rightpart(tl + h * 13 / 20, xl + h * (k1  * (246121993. / 1340847787.) + k4 * ((-37695042795.) / 15268766246.) + k5 * ((-309121744.) / 1061227803.) +
                k6 * ((-12992083.) / 490766935.) + k7 * (6005943493. / 2108947869.) + k8 * (393006217. / 1396673457) + k9 * (123872331. / 1001029789.) ), eps);

        k11 = rightpart(tl + h * (1201146811. / 1299019798.), xl + h * (k1  * ((-1028468189.) / 846180014.) + k4 * (8478235783. / 508512852.) + k5 * (1311729495. / 1432422823.) +
                k6 * ((-10304129995.) / 1701304382.) + k7 * ((-48777925059.) / 3047939560.) + k8 * (15336726248. / 1032824649.) + k9 * ((-45442868181.) / 3398467696.) + 
                k10 * (3065993473. / 597172653.) ), eps);

        k12 = rightpart(tl + h, xl + h * (k1  * (185892177. / 718116043.) + k4 * ((-3185094517.) / 667107341.) + k5 * ((-477755414.) / 1098053517.) + 
                 k6 * ((-703635378.) / 230739211.) + k7 * (5731566787. / 1027545527.) + k8 * (5232866602. / 850066563.) + k9 * ((-4093664535.) / 808688257.) + 
                k10 * (3962137247. / 1805957418.) + k11 * (65686358. / 487910083.) ), eps);

        k13 = rightpart(tl + h, xl + h * (k1  * (403863854. / 491063109.) + k4 * ((-5068492393.) / 434740067.) + k5 * ((-411421997.) / 543043805.) + 
                 k6 * (652783627. / 914296604.) + k7 * (11173962825. / 925320556.) + k8 * ((-13158990841.) / 6184727034.) + k9 * (3936647629. / 1978049680.) + 
                k10 * ((-160528059.) / 685178525.) + k11 * (248638103. / 1413531060.) ), eps);

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
            chi_2_vec.append(xl[6]);

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
                chi_2_vec.append(xl[6]);

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
}
