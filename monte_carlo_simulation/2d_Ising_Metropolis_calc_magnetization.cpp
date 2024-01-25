#include <iostream>
#include <cmath>
#include <fstream>
const long int niter = 100000;
const int L = 128;
const int nx = L; // number of sites along x-direction
const int ny = L; // number of sites along y-direction
const int monte_calro_step = niter * nx * ny;
const double coupling_J = 1.0;
const int nconf = 80;
const double t_start = 1.8;
const int nskip = nx * ny * 10; // Frequency of measurement
const int nconfig = 1;

double calc_action_change(const int spin[nx][ny], const double coupling_J, const double temperature, const int ix, const int iy)
{
    double action_change = 0e0;

    int ixp1 = (ix + 1) % nx;
    int iyp1 = (iy + 1) % ny;
    int ixm1 = (ix - 1 + nx) % nx;
    int iym1 = (iy - 1 + ny) % ny;

    double sum_change = 2 * spin[ix][iy] * spin[ixp1][iy] + 2 * spin[ix][iy] * spin[ix][iyp1] + 2 * spin[ix][iy] * spin[ixm1][iy] + 2 * spin[ix][iy] * spin[ix][iym1];

    action_change = sum_change * coupling_J / temperature;

    return action_change;
}

double calc_total_spin(const int spin[nx][ny])
{
    double total_spin = 0;
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            total_spin += spin[ix][iy];
        }
    }
    if (total_spin < 0)
    {
        total_spin = -total_spin;
    }
    return total_spin;
}

double calc_magnetization(const int total_spin, const int total_sampling)
{
    double magnetization = 0;
    magnetization = (double)total_spin / (double)total_sampling;
    if (magnetization < 0)
    {
        magnetization = -magnetization;
    }
    return magnetization;
}

int main()
{
    double temperature[nconf + 1];
    double sum = t_start;
    for (int i = 0; i < nconf + 1; i++)
    {
        temperature[i] = sum;
        sum += 0.01;
    }
    std::ofstream outputfile("output/2d_Ising_L" + std::to_string(L) + "_magnetization_metropolis.txt");
    for (int conf = 0; conf < nconf + 1; conf++)
    {
        double T = temperature[conf];
        // 初期化
        int count = 0;
        int total_spin_sum = 0;
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
            std::ifstream inputconfig("output/2d_Ising_Metropolis_output_config.txt");
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
        // 各温度でのモンテカルロシミュレーション
        for (long int iter = 0; iter != monte_calro_step; iter++)
        {
            double rand_site = (double)rand() / RAND_MAX;
            rand_site = rand_site * nx * ny;
            int ix = (int)rand_site / ny;
            int iy = (int)rand_site % ny;
            double metropolis = (double)rand() / RAND_MAX;
            double action_change = calc_action_change(spin, coupling_J, T, ix, iy);
            if (exp(-action_change) > metropolis)
            {
                // accept
                spin[ix][iy] = -spin[ix][iy]; // flip
            }
            else
            {
                // reject
            }
            if (iter > 1000 * nx * ny && (iter + 1) % nskip == 0)
            {
                total_spin_sum += calc_total_spin(spin);
                count++;
            }
        }
        std::cout << count << std::endl;
        double magnetization = calc_magnetization(total_spin_sum, count);
        double m = magnetization / (nx * ny);
        std::cout << std::fixed << std::setprecision(4)
                  << T << "   "
                  << m
                  << std::endl;
        outputfile << std::fixed << std::setprecision(4)
                   << T << "   "
                   << m
                   << std::endl;
    }
    outputfile.close();
    return 0;
}