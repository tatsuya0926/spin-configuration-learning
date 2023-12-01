#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

const int monte_carlo_step = 2000;
const int niter = 1000000;
const double step_size = 0.1e0;
const double cooling_rate = 0.995;

double calc_f(double x)
{
    double fx = (x - 1e0) * (x - 1e0) * ((x + 1e0) * (x + 1e0) + 0.01e0);

    return fx;
}

int main(void)
{
    int naccept[monte_carlo_step];
    double x[monte_carlo_step];
    double temperature = 5e0;

    srand((unsigned)time(NULL));

    /*********************************/
    /*********** 初期状態 ************/
    /*********************************/
    for (int i = 0; i < monte_carlo_step; i++)
    {
        x[i] = 0e0;
        naccept[i] = 0;
    }
    /**************/
    /** Main処理 **/
    /**************/
    for (int iter = 1; iter < niter + 1; iter++)
    {
        for (int step = 0; step < monte_carlo_step; step++)
        {
            double backup_x = x[step];
            double current_obj = calc_f(x[step]);

            // double dx = (double)rand() / RAND_MAX;
            // dx = (dx - 0.5e0) * step_size * 2e0;
            // x[step] = x[step] + dx;

            // ボックス・ミュラー変換
            double p = (double)rand() / RAND_MAX;
            double q = (double)rand() / RAND_MAX;
            x[step] = x[step] + temperature * (sqrt(-2 * log(p)) * cos(2 * M_PI * q));

            double new_obj = calc_f(x[step]);
            double metropolis = (double)rand() / RAND_MAX;
            if (exp((current_obj - new_obj) / temperature) > metropolis)
                /* accept */
                naccept[step] = naccept[step] + 1;
            else
                /* reject */
                x[step] = backup_x;
        }
        temperature = cooling_rate * temperature;
        /***************/
        /* data output */
        /***************/
        printf("%f    %f    %f\n", x[19], x[199], x[1999]);
    }
}
