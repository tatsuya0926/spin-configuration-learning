#include <iostream>
#include <cmath>
#include <fstream>
const long int monte_carlo_step = 100000;
const int L = 64;
const int nx = L; // number of sites along x-direction
const int ny = L; // number of sites along y-direction
const int niter = monte_carlo_step * nx * ny;
const int Q = 5;
const double coupling_J = 1.0;
const double temperature = 5.0;
const int nconfig = 1;
const int nskip = nx * ny * 100; // Frequency of measurement

double kronecker_delta(const int spin_1, const int spin_2)
{
    if (spin_1 == spin_2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

double calc_action_change(const int spin[nx][ny], const int next_spin, const double coupling_J, const double temperature, const int ix, const int iy)
{
    double action_change = 0e0;

    int ixp1 = (ix + 1) % nx;      // ixp1=ix+1; be careful about the boundary condition.
    int iyp1 = (iy + 1) % ny;      // iyp1=iy+1; be careful about the boundary condition.
    int ixm1 = (ix - 1 + nx) % nx; // ixm1=ix-1; be careful about the boundary condition.
    int iym1 = (iy - 1 + ny) % ny; // iym1=iy-1; be careful about the boundary condition.

    double sum_change = kronecker_delta(spin[ix][iy], spin[ixp1][iy]) +
                        kronecker_delta(spin[ix][iy], spin[ix][iyp1]) +
                        kronecker_delta(spin[ix][iy], spin[ixm1][iy]) +
                        kronecker_delta(spin[ix][iy], spin[ix][iym1]) -
                        kronecker_delta(next_spin, spin[ixp1][iy]) -
                        kronecker_delta(next_spin, spin[ix][iyp1]) -
                        kronecker_delta(next_spin, spin[ixm1][iy]) -
                        kronecker_delta(next_spin, spin[ix][iym1]);

    action_change = sum_change * coupling_J / temperature;

    return action_change;
}

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
        double action_change = calc_action_change(spin, coupling_J, next_spin, temperature, ix, iy);
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
    std::ofstream outputconfig("../output/2d_Potts_q=" + std::to_string(Q) + "_output_config.txt");
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
