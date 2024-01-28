#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
const long int monte_carlo_step = 100000;
const int L = 64;
const int nx = L; // number of sites along x-direction
const int ny = L; // number of sites along y-direction
const int niter = monte_carlo_step * nx * ny;
const double coupling_J = 1.0;
const int nconf = 30;
const int ndata = 1000;
const double t_start = 2.1;
const int nskip = nx * ny * 100; // Frequency of measurement
const int nconfig = 0;

double calc_action_change(const int spin[nx][ny], const double coupling_J, const double temperature, const int ix, const int iy)
{
    double action_change = 0.0;
    // double sum_change=0;

    int ixp1 = (ix + 1) % nx;      // ixp1=ix+1; be careful about the boundary condition.
    int iyp1 = (iy + 1) % ny;      // iyp1=iy+1; be careful about the boundary condition.
    int ixm1 = (ix - 1 + nx) % nx; // ixm1=ix-1; be careful about the boundary condition.
    int iym1 = (iy - 1 + ny) % ny; // iym1=iy-1; be careful about the boundary condition.

    double sum_change = 2 * spin[ix][iy] * spin[ixp1][iy] +
                        2 * spin[ix][iy] * spin[ix][iyp1] +
                        2 * spin[ix][iy] * spin[ixm1][iy] +
                        2 * spin[ix][iy] * spin[ix][iym1];

    action_change = sum_change * coupling_J / temperature;

    return action_change;
}

int main()
{
    double temperature[nconf + 1];
    double sum = t_start;
    for (int i = 0; i < nconf + 1; i++)
    {
        temperature[i] = sum;
        sum += 0.01;
        std::cout << temperature[i] << std::endl;
    }
    for (int conf = 0; conf < nconf + 1; conf++)
    {
        double T = temperature[conf];
        int data_num = 0;
        int spin[nx][ny];
        srand((unsigned)time(NULL));
        // 初期化
        if (nconfig == 0)
        {
            std::ifstream inputconfig("input/2d_Ising_output_config.txt");
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
        for (long int iter = 0; iter != niter; iter++)
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
            if (iter > 1000 * nx * ny && (iter + 1) % nskip == 0 && data_num < ndata)
            {
                std::ofstream outputfile("../txtfile/2d_Ising/L" + std::to_string(L) + "T" + std::to_string(conf) + "_" + std::to_string(data_num+ndata) + ".txt");
                for (int ix = 0; ix != nx; ix++)
                {
                    for (int iy = 0; iy != ny; iy++)
                    {
                        outputfile << ix << ' ' << iy << ' ' << spin[ix][iy] << ' ' << std::endl;
                    }
                }
                outputfile.close();
                data_num++;
            }
        }
    }
    return 0;
}