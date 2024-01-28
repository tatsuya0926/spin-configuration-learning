#include <iostream>
#include <cmath>
#include <fstream>
const double pi = 3.141592653589793;
const long int monte_carlo_step = 100000;
const int L = 64;
const int nx = L; // number of sites along x-direction
const int ny = L; // number of sites along y-direction
const int niter = monte_carlo_step * nx * ny;
const int Q = 6;
const double coupling_J = 1.0;
const double temperature = 5.0;
const int nconfig = 1;
const int nskip = nx * ny * 100; // Frequency of measurement

double calc_action_change(const int spin[nx][ny], const int next_spin, const double coupling_J, const double temperature, const int ix, const int iy)
{

    double action_change = 0e0;

    int ixp1 = (ix + 1) % nx;      // ixp1=ix+1; be careful about the boundary condition.
    int iyp1 = (iy + 1) % ny;      // iyp1=iy+1; be careful about the boundary condition.
    int ixm1 = (ix - 1 + nx) % nx; // ixm1=ix-1; be careful about the boundary condition.
    int iym1 = (iy - 1 + ny) % ny; // iym1=iy-1; be careful about the boundary condition.

    int sum_change = std::cos(2 * pi * spin[ix][iy] / Q - 2 * pi * spin[ixp1][iy] / Q) +
                     std::cos(2 * pi * spin[ix][iy] / Q - 2 * pi * spin[ix][iyp1] / Q) +
                     std::cos(2 * pi * spin[ix][iy] / Q - 2 * pi * spin[ixm1][iy] / Q) +
                     std::cos(2 * pi * spin[ix][iy] / Q - 2 * pi * spin[ix][iym1] / Q) -
                     std::cos(2 * pi * next_spin / Q - 2 * pi * spin[ixp1][iy] / Q) -
                     std::cos(2 * pi * next_spin / Q - 2 * pi * spin[ix][iyp1] / Q) -
                     std::cos(2 * pi * next_spin / Q - 2 * pi * spin[ixm1][iy] / Q) -
                     std::cos(2 * pi * next_spin / Q - 2 * pi * spin[ix][iym1] / Q);

    action_change = sum_change * coupling_J / temperature;

    return action_change;
}
/*************************************/
/*** Calculation of the total spin ***/
/*************************************/
int calc_total_spin(const int spin[nx][ny])
{

    int total_spin = 0;

    for (int ix = 0; ix != nx; ix++)
    {
        for (int iy = 0; iy != ny; iy++)
        {
            total_spin = total_spin + spin[ix][iy];
        }
    }

    return total_spin;
}
/************/
/*** Main ***/
/************/
int main()
{
    int spin[nx][ny];
    srand((unsigned)time(NULL));
    /*********************************/
    /* Set the initial configuration */
    /*********************************/
    for (int ix = 0; ix != nx; ix++)
    {
        for (int iy = 0; iy != ny; iy++)
        {
            spin[ix][iy] = nconfig;
        }
    }
    // std::ofstream outputfile("output.txt");
    int naccept = 0; // counter for the number of acceptance

    for (long int iter = 0; iter != niter; iter++)
    {
        // choose a point randomly.
        double rand_site = (double)rand() / RAND_MAX;
        rand_site = rand_site * nx * ny;
        int ix = (int)rand_site / ny;
        int iy = (int)rand_site % ny;
        double metropolis = (double)rand() / RAND_MAX;
        int next_spin = rand() % Q;
        double action_change = calc_action_change(spin, next_spin, coupling_J, temperature, ix, iy);
        if (exp(-action_change) > metropolis)
        {
            // accept
            spin[ix][iy] = next_spin;
            naccept = naccept + 1;
        }
        else
        {
            // reject
        }
        // int total_spin = calc_total_spin(spin);
        // double energy = calc_action(spin, coupling_J, coupling_h, temperature) * temperature;
        // if ((iter + 1) % nskip == 0)
        // {
        //     std::cout << std::fixed << std::setprecision(4)
        //               << total_spin << "   "
        //               << energy << "   "
        //               << (double)naccept / (iter + 1) << std::endl;
        //     outputfile << std::fixed << std::setprecision(4)
        //                << total_spin << "   "
        //                << energy << "   "
        //                << (double)naccept / (iter + 1) << std::endl;
        // }
    }
    // outputfile.close();
    std::ofstream outputconfig("../output/2d_Clock_q=" + std::to_string(Q) + "_output_config.txt");
    for (int ix = 0; ix != nx; ix++)
    {
        for (int iy = 0; iy != ny; iy++)
        {
            outputconfig << ix << ' ' << iy << ' ' << spin[ix][iy]
                         << ' ' << std::endl;
        }
    }
    outputconfig.close();
    return 0;
}
