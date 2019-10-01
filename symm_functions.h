#ifndef FUNCTIONS_H
#define FUNCTIONS_h

#include "include.h"

double theta            (double t, double D3, double W3, double D6);
double beta             (double t1, double t, double D3, double W3, double D6);
double gamma_W_c        (double t1, double t, double D3, double W3, double D6);
double gamma_W_s        (double t1, double t, double D3, double W3, double D6);
double Gamma_D_1        (double t1, double t, double D3, double W3, double D6);
double Gamma_D_2        (double t1, double t, double D3, double W3, double D6);
double Gamma_W_1        (double t1, double t, double D3, double W3, double D6);
double Gamma_W_2        (double t1, double t, double D3, double W3, double D6);
double gamma_s_plus_c   (double t1, double t, double D3, double W3, double D6);
double gamma_s_minus_c  (double t1, double t, double D3, double W3, double D6);

double u_control_1      (double W_control_1, double W_control_2, double W_control_3);
double u_control_2      (double W_control_1, double W_control_2, double W_control_3);
double u_control_3      (double W_control_1, double W_control_2, double W_control_3);
Vector<3> W_control     (Vector<3> u_control);

#endif //FUNCTIONS_H
