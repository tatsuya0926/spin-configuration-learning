#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

const int nbeta = 2000;
const int niter = 1000000;
const double step_size = 0.1e0;
const double dbeta = 0.5e0;

double calc_f(double x)
{
    double fx = (x - 1e0) * (x - 1e0) * ((x + 1e0) * (x + 1e0) + 0.01e0);

    return fx;
}

int main(void)
{
    int naccept[nbeta];
    double x[nbeta];
    double beta[nbeta];

    srand((unsigned)time(NULL));

    /*********************************/
    /*********** 初期状態 ************/
    /*********************************/
    for (int ibeta = 0; ibeta < nbeta; ibeta++)
    {
        x[ibeta] = 0e0;
        beta[ibeta] = (double)(ibeta + 1) * dbeta;
        naccept[ibeta] = 0;
    }
    /**************/
    /** Main処理 **/
    /**************/
    for (int iter = 1; iter < niter + 1; iter++)
    {
        /* 各レプリカにおけるモンテカルロ */
        for (int ibeta = 0; ibeta < nbeta; ibeta++)
        {
            double backup_x = x[ibeta];
            double current_obj = calc_f(x[ibeta]) * beta[ibeta];

            double dx = (double)rand() / RAND_MAX;
            dx = (dx - 0.5e0) * step_size * 2e0;
            x[ibeta] = x[ibeta] + dx;

            double new_obj = calc_f(x[ibeta]) * beta[ibeta];
            double metropolis = (double)rand() / RAND_MAX;
            if (exp(current_obj - new_obj) > metropolis)
                /* accept */
                naccept[ibeta] = naccept[ibeta] + 1;
            else
                /* reject */
                x[ibeta] = backup_x;
        }
        /* レプリカ交換処理 */
        for (int ibeta = 0; ibeta < nbeta - 1; ibeta++)
        {
            double current_obj = calc_f(x[ibeta]) * beta[ibeta] + calc_f(x[ibeta + 1]) * beta[ibeta + 1];
            double new_obj = calc_f(x[ibeta]) * beta[ibeta + 1] + calc_f(x[ibeta + 1]) * beta[ibeta];
            double metropolis = (double)rand() / RAND_MAX;
            if (exp(current_obj - new_obj) > metropolis)
            {
                double backup_x = x[ibeta];
                x[ibeta] = x[ibeta + 1];
                x[ibeta + 1] = backup_x;
            }
        }

        /***************/
        /* data output */
        /***************/
        printf("%f    %f    %f\n", x[19], x[199], x[1999]);
    }
}
