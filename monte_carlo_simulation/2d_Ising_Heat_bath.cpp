#include <iostream>
#include <cmath>
#include <fstream>
const long int niter = 1000000;
const int nx = 64; // number of sites along x-direction
const int ny = 64; // number of sites along y-direction
const double coupling_J = 1.0;
const double coupling_h = 0.1;
const double temperature = 5.0;
const int nskip = 100; // Frequency of measurement
const int nconfig = 1; // 0 -> read 'input_config.txt'; 1 -> all up; -1 -> all down

double calc_action(const int spin[nx][ny], const double coupling_J, const double coupling_h, const double temperature)
{
    double action = 0e0;
    int sum1 = 0;
    int sum2 = 0;
    for (int ix = 0; ix != nx; ix++)
    {
        int ixp1 = (ix + 1) % nx;
        for (int iy = 0; iy != ny; iy++)
        {
            int iyp1 = (iy + 1) % ny;
            sum1 = sum1 + spin[ix][iy];
            sum2 = sum2 + spin[ix][iy] * spin[ixp1][iy] + spin[ix][iy] * spin[ix][iyp1];
        }
    }
    action = (sum2 * coupling_J + sum1 * coupling_h) / temperature * (-1e0);

    return action;
}

double heat_bath_probability(const int spin[nx][ny], const double coupling_J, const double coupling_h, const double temperature, const int ix, const int iy)
{
    double action_change = 0e0;
    double temp = 0e0;
    double Ep, Em; // E_+, E_-

    int ixp1 = (ix + 1) % nx;
    int iyp1 = (iy + 1) % ny;
    int ixm1 = (ix - 1 + nx) % nx;
    int iym1 = (iy - 1 + ny) % ny;

    temp = coupling_h;
    temp = temp + (spin[ixp1][iy] + spin[ix][iyp1] + spin[ixm1][iy] + spin[ix][iym1]) * coupling_J;

    Ep = temp / temperature * (-1e0);
    Em = temp / temperature;

    double ratio = exp(-Ep) / (exp(-Ep) + exp(-Em));

    return ratio;
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

int calc_total_plus_spin(const int spin[nx][ny])
{
    int total_spin = 0;
    for (int ix = 0; ix != nx; ix++)
    {
        for (int iy = 0; iy != ny; iy++)
        {
            if (spin[ix][iy] == 1)
            {
                total_spin = total_spin + 1;
            }
        }
    }
    return total_spin;
}

int main()
{
    int spin[nx][ny];
    srand((unsigned)time(NULL));
    if (nconfig == 1)
    {
        for (int ix = 0; ix != nx; ix++)
        {
            for (int iy = 0; iy != ny; iy++)
            {
                spin[ix][iy] = 1;
            }
        }
    }
    if (nconfig == -1)
    {
        for (int ix = 0; ix != nx; ix++)
        {
            for (int iy = 0; iy != ny; iy++)
            {
                spin[ix][iy] = -1;
            }
        }
    }
    if (nconfig == 0)
    {
        std::ifstream inputconfig("input_config.txt");
        if (!inputconfig)
        {
            std::cout << "inputfile not found" << std::endl;
            exit(1);
        }
        for (int ix = 0; ix != nx; ix++)
        {
            for (int iy = 0; iy != ny; iy++)
            {
                inputconfig >> ix >> iy >> spin[ix][iy];
            }
        }
        inputconfig.close();
    }
    std::ofstream outputfile("output/Heat_bath/2d_Ising_Heat_bath_output_t5.txt");
    int count = 0;

    for (long int iter = 0; iter != niter; iter++)
    {
        double rand_site = (double)rand() / RAND_MAX;
        rand_site = rand_site * nx * ny;
        int ix = (int)rand_site / ny;
        int iy = (int)rand_site % ny;
        double metropolis = (double)rand() / RAND_MAX;

        double ratio = heat_bath_probability(spin, coupling_J, coupling_h, temperature, ix, iy);
        if (ratio > metropolis)
        {
            spin[ix][iy] = 1;
        }
        else
        {
            spin[ix][iy] = -1;
        }

        int total_spin = calc_total_spin(spin);
        int total_plus_spin = calc_total_plus_spin(spin);
        double energy = calc_action(spin, coupling_J, coupling_h, temperature) * temperature;
        // if (iter % nskip == 0)
        // {
        //     std::ofstream outputconfig("output/fig/2d_Ising_Heat_bath_output_config_" + std::to_string(count) + ".txt");
        //     for (int ix = 0; ix != nx; ix++)
        //     {
        //         for (int iy = 0; iy != ny; iy++)
        //         {
        //             outputconfig << ix << ' ' << iy << ' ' << spin[ix][iy] << ' ' << std::endl;
        //         }
        //     }
        //     count++;
        //     outputconfig.close();
        // }

        if ((iter + 1) % nskip == 0)
        {
            std::cout << std::fixed << std::setprecision(4)
                      << count * nskip << "   "
                      << total_spin << "   "
                      << total_plus_spin << "   "
                      << energy << "   " << std::endl;
            outputfile << std::fixed << std::setprecision(4)
                       << count * nskip << "   "
                       << total_spin << "   "
                       << total_plus_spin << "   "
                       << energy << "   " << std::endl;
            count++;
        }
    }
    outputfile.close();
    return 0;
}