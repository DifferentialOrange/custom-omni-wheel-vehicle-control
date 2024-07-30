#include "vector.h"
#include "parameters.h"
#include "include.h"
#include <QDebug>

Vector<6> custom_control_find(Vector<6> control, double t_sw, double T, Vector<6> initial_values,
                              Vector<6> final_values)
{
    auto alpha = parameters::symmetrical::alpha;
    auto beta = parameters::symmetrical::beta;
    auto delta = parameters::symmetrical::delta;
    double Delta = parameters::symmetrical::Delta;
    double Lambda = parameters::symmetrical::Lambda;

    double angle_step = precision::angle_step;
    double N_step = precision::N_step;

    for (unsigned i = 0; i < 3; ++i) {
        qDebug() << "alpha" << i + 1 << " STARTING\n";
        while(fabs(parameters::final::alpha[i] - alpha[i]) > precision::double_EPS)
        {
            qDebug() << "alpha" << i + 1 << " COMPUTING to " << parameters::final::alpha[i] << "\n";
            auto pot_alpha = alpha;
            auto pot_beta = beta;
            auto pot_delta = delta;
            auto pot_Delta = Delta;
            auto pot_Lambda = Lambda;

            if (pot_alpha[i] > parameters::final::alpha[i])
            {
                auto step = std::min(angle_step,
                                     fabs(parameters::final::alpha[i] - pot_alpha[i]));
                pot_alpha[i] -= step;
                pot_beta[i] -= step;

                if (i == 0)
                {
                    pot_alpha[2] += step;
                    pot_beta[2] += step;
                }
            }
            else
            {
                auto step = std::min(angle_step,
                                     fabs(parameters::final::alpha[i] - pot_alpha[i]));
                pot_alpha[i] += step;
                pot_beta[i] += step;

                if (i == 0)
                {
                    pot_alpha[2] -= step;
                    pot_beta[2] -= step;
                }
            }

            auto ans = control_solve(0, T, initial_values, final_values, control,
                                t_sw, pot_alpha, pot_beta, pot_delta, pot_Delta, pot_Lambda);

            if (ans.first)
            {
                alpha = pot_alpha;
                beta = pot_beta;
                delta = pot_delta;
                Delta = pot_Delta;
                Lambda = pot_Lambda;

                control = ans.second;

                qDebug() << "alpha" << i + 1 << " " << alpha[i] << " OK\n";
            }
            else
            {
                qDebug() << "angle_step " << angle_step << "\n";
                angle_step /= 2.0;
            }
        }
    }


    for (unsigned i = 0; i < 3; ++i) {
        qDebug() << "beta" << i + 1 << " STARTING\n";
        while(fabs(parameters::final::beta[i] - beta[i]) > precision::double_EPS)
        {
            qDebug() << "beta" << i + 1 << " COMPUTING to " << parameters::final::beta[i] << "\n";
            auto pot_alpha = alpha;
            auto pot_beta = beta;
            auto pot_delta = delta;
            auto pot_Delta = Delta;
            auto pot_Lambda = Lambda;

            if (pot_beta[i] > parameters::final::beta[i])
            {
                auto step = std::min(angle_step,
                                     fabs(parameters::final::beta[i] - pot_beta[i]));
                pot_beta[i] -= step;

                if (i == 0)
                {
                    pot_beta[2] += step;
                }
            }
            else
            {
                auto step = std::min(angle_step,
                                     fabs(parameters::final::beta[i] - pot_beta[i]));
                pot_beta[i] += step;

                if (i == 0)
                {
                    pot_beta[2] -= step;
                }
            }

            auto ans = control_solve(0, T, initial_values, final_values, control,
                                t_sw, pot_alpha, pot_beta, pot_delta, pot_Delta, pot_Lambda);

            if (ans.first)
            {
                alpha = pot_alpha;
                beta = pot_beta;
                delta = pot_delta;
                Delta = pot_Delta;
                Lambda = pot_Lambda;

                control = ans.second;

                qDebug() << "beta" << i + 1 << " " << beta[i] << " OK\n";
            }
            else
            {
                qDebug() << "angle_step " << angle_step << "\n";
                angle_step /= 2.0;
            }
        }
    }

    for (unsigned i = 0; i < 3; ++i) {
        qDebug() << "delta" << i + 1 << " STARTING\n";
        while(fabs(parameters::final::delta[i] - delta[i]) > precision::double_EPS)
        {
            qDebug() << "delta" << i + 1 << " COMPUTING to " << parameters::final::delta[i] << "\n";
            auto pot_alpha = alpha;
            auto pot_beta = beta;
            auto pot_delta = delta;
            auto pot_Delta = Delta;
            auto pot_Lambda = Lambda;

            if (pot_delta[i] > parameters::final::delta[i])
            {
                auto step = std::min(N_step,
                                     fabs(parameters::final::delta[i] - pot_delta[i]));
                pot_delta[i] -= step;

                if (i == 0)
                {
                    pot_delta[2] -= step;
                }

                if (i == 1)
                {
                    pot_Delta -= step;
                }
            }
            else
            {
                auto step = std::min(N_step,
                                     fabs(parameters::final::delta[i] - pot_delta[i]));
                pot_delta[i] += step;

                if (i == 0)
                {
                    pot_delta[2] += step;
                }

                if (i == 1)
                {
                    pot_Delta += step;
                }
            }

            auto ans = control_solve(0, T, initial_values, final_values, control,
                                t_sw, pot_alpha, pot_beta, pot_delta, pot_Delta, pot_Lambda);

            if (ans.first)
            {
                alpha = pot_alpha;
                beta = pot_beta;
                delta = pot_delta;
                Delta = pot_Delta;
                Lambda = pot_Lambda;

                control = ans.second;

                qDebug() << "delta" << i + 1 << " " << delta[i] << "OK\n";
            }
            else
            {
                N_step /= 2.0;
                qDebug() << "N_step " << N_step << "\n";
            }
        }
    }

    qDebug() << "Delta STARTING\n";
    while(fabs(parameters::final::Delta - Delta) > precision::double_EPS)
    {
        qDebug() << "Delta COMPUTING to " << parameters::final::Delta << "\n";
        auto pot_alpha = alpha;
        auto pot_beta = beta;
        auto pot_delta = delta;
        auto pot_Delta = Delta;
        auto pot_Lambda = Lambda;

        if (pot_Delta > parameters::final::Delta)
        {
            pot_Delta -= std::min(N_step,
                                  fabs(parameters::final::Delta - pot_Delta));
        }
        else
        {
            pot_Delta += std::min(N_step,
                                  fabs(parameters::final::Delta - pot_Delta));
        }

        auto ans = control_solve(0, T, initial_values, final_values, control,
                            t_sw, pot_alpha, pot_beta, pot_delta, pot_Delta, pot_Lambda);

        if (ans.first)
        {
            alpha = pot_alpha;
            beta = pot_beta;
            delta = pot_delta;
            Delta = pot_Delta;
            Lambda = pot_Lambda;

            control = ans.second;
        }
        else
        {
            N_step /= 2.0;
            qDebug() << "N_step " << N_step << "\n";
        }

        qDebug() << "Delta " << Delta << " OK\n";
    }

    qDebug() << "Lambda STARTING\n";
    while(fabs(parameters::final::Lambda - Lambda) > precision::double_EPS)
    {
        qDebug() << "Lambda COMPUTING to " << parameters::final::Lambda << "\n";
        auto pot_alpha = alpha;
        auto pot_beta = beta;
        auto pot_delta = delta;
        auto pot_Delta = Delta;
        auto pot_Lambda = Lambda;

        if (pot_Lambda > parameters::final::Lambda)
        {
            pot_Lambda -= std::min(N_step,
                                  fabs(parameters::final::Lambda - pot_Lambda));
        }
        else
        {
            pot_Lambda += std::min(N_step,
                                  fabs(parameters::final::Lambda - pot_Lambda));
        }

        auto ans = control_solve(0, T, initial_values, final_values, control,
                            t_sw, pot_alpha, pot_beta, pot_delta, pot_Delta, pot_Lambda);

        if (ans.first)
        {
            alpha = pot_alpha;
            beta = pot_beta;
            delta = pot_delta;
            Delta = pot_Delta;
            Lambda = pot_Lambda;

            control = ans.second;
        }
        else
        {
            N_step /= 2.0;
            qDebug() << "N_step " << N_step << "\n";
        }

        qDebug() << "Lambda " << Lambda << " OK\n";
    }


    return control;
}
