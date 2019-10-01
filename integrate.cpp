#include <cmath>
#include <functional>

static const double EPS = 1e-9;

double DOPRI8_par_integrate(double t_left, double t_right, std::function<double(double, double, double, double, double)> int_function,
                                 double D3, double W3, double D6)
{
    if (t_right - t_left < 1e-15)
        return 0;

    double h = (t_right - t_left) / 1e7;
    double h_new;

    double tl = t_left;
    double xl = 0; 

    double stepx; 

    double errx; 
    double err, coef;
    double coefmax = 5;
    bool last_flag = false;

    double k1, /*k2, k3, k4, k5,*/ k6, k7, k8, k9, k10, k11, k12, k13;


    while (tl + h < t_right || last_flag)
    {

        k1 = int_function(t_left, tl, D3, W3, D6);

        /*k2 = int_function(t_left, tl + h / 18, D3, W3, D6);

        k3 = int_function(t_left, tl + h / 12, D3, W3, D6);

        k4 = int_function(t_left, tl + h / 8, D3, W3, D6);

        k5 = int_function(t_left, tl + h * 5 / 16, D3, W3, D6);*/

        k6 = int_function(t_left, tl + h * 3 / 8, D3, W3, D6);

        k7 = int_function(t_left, tl + h * 59. / 400., D3, W3, D6);

        k8 = int_function(t_left, tl + h * 93. / 200., D3, W3, D6);

        k9 = int_function(t_left, tl + h * (5490023248. / 9719169821.), D3, W3, D6);

        k10 = int_function(t_left, tl + h * 13 / 20, D3, W3, D6);

        k11 = int_function(t_left, tl + h * (1201146811. / 1299019798.), D3, W3, D6);

        k12 = int_function(t_left, tl + h, D3, W3, D6);

        k13 = int_function(t_left, tl + h, D3, W3, D6);


        stepx = h * (k1 * (14005451. / 335480064.) + k6 * ((-59238493.) / 1068277825.)
                     + k7 * (181606767. / 758867731.) + k8 * (561292985. / 797845732.) +
                     k9 * ((-1041891430.) / 1371343529.) + k10 * (760417239. / 1151165299.)
                     + k11 * (118820643. / 751138087.) + k12 * ((-528747749.) / 2220607170.) + k13 / 4 );


        errx = h * (k1 * (13451932. / 455176623.) + k6 * ((-808719846.) / 976000145.)
                    + k7 * (1757004468. / 5645159321.) + k8 * (656045339. / 265891186.) +
                    k9 * ((-3867574721.) / 1518517206.) + k10 * (465885868. / 322736535.)
                    + k11 * (53011238. / 667516719.) + k12 * 2 / 45 );

        //step correction

        err = fabs(stepx - errx);

        if (err < 1e-15)
        {
            tl += h;
            xl += stepx;

            coefmax = 5;

            h *= 2;
        }
        else
        {
            coef = pow(EPS / err, 1.0 / (8 + 1));

            if (coef > coefmax)
                coef = coefmax;
            else if (coef < 0.2)
                coef = 0.2;

            h_new = 0.9 * h * coef;
    
            if (EPS > err)
            {            
                tl += h;
                xl += stepx;
    
                coefmax = 5;
            }
            else
                coefmax = 1;

            h = h_new;
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
