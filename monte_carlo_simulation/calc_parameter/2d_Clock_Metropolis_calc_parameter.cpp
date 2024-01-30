#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
const double pi = 3.141592653589793;
const long int niter = 100000;
const int L = 64;
const int nx = L; // number of sites along x-direction
const int ny = L; // number of sites along y-direction
const int monte_carlo_step = niter * nx * ny;
const int Q = 4;
const double coupling_J = 1.0;
const int nconf = 80;
const double t_start = 1.8;
const int nskip = nx * ny * 10; // Frequency of measurement
const int nconfig = 1;

double calc_energy(const int spin[nx][ny], const double coupling_J, const double temperature)
{
    double action = 0.0;
    double sum = 0;
    for (int ix = 0; ix != nx; ix++)
    {
        int ixp1 = (ix + 1) % nx;
        for (int iy = 0; iy != ny; iy++)
        {
            int iyp1 = (iy + 1) % ny;
            sum = sum + std::cos(2 * pi * spin[ix][iy] / Q - 2 * pi * spin[ixp1][iy] / Q) + std::cos(2 * pi * spin[ix][iy] / Q - 2 * pi * spin[ix][iyp1] / Q);
        }
    }
    action = sum * coupling_J / temperature * (-1e0);

    return action;
}

double calc_action_change(const int spin[nx][ny], const int next_spin, const double coupling_J, const double temperature, const int ix, const int iy)
{
    double action_change = 0e0;
    // double sum_change=0;

    int ixp1 = (ix + 1) % nx;      // ixp1=ix+1;
    int iyp1 = (iy + 1) % ny;      // iyp1=iy+1;
    int ixm1 = (ix - 1 + nx) % nx; // ixm1=ix-1;
    int iym1 = (iy - 1 + ny) % ny; // iym1=iy-1;

    double sum_change = std::cos(2 * pi * spin[ix][iy] / Q - 2 * pi * spin[ixp1][iy] / Q) +
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

double calc_squared_magnetization(const int spin[nx][ny])
{
    std::vector<double> m(Q, 0.0);
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            m[spin[ix][iy]] += 1;
        }
    }
    for (int i = 0; i < Q; i++)
    {
        m[i] /= static_cast<double>(nx * ny);
    }
    double m2 = 0.0;
    for (int i = 0; i < Q; i++)
    {
        m2 += m[i] * m[i];
    }
    for (int i = 0; i < Q - 1; i++)
    {
        for (int j = i + 1; j < Q; j++)
        {
            m2 -= 2.0 * m[i] * m[j] / (Q - 1);
        }
    }
    return m2;
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
    std::ofstream outputfile("output/2d_Clock_L" + std::to_string(L) + "_q=" + std::to_string(Q) + "_parameter_metropolis.txt");
    for (int conf = 0; conf < nconf + 1; conf++)
    {
        double T = temperature[conf];
        // 初期化
        int count = 0;
        int total_m2 = 0;
        int total_m4 = 0;
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
            std::ifstream inputconfig("output/2d_Clock_Metropolis_output_config.txt");
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
        for (long int iter = 0; iter != monte_carlo_step; iter++)
        {
            double rand_site = (double)rand() / RAND_MAX;
            rand_site = rand_site * nx * ny;
            int ix = (int)rand_site / ny;
            int iy = (int)rand_site % ny;
            double metropolis = (double)rand() / RAND_MAX;
            int next_spin = rand() % Q;
            double action_change = calc_action_change(spin, next_spin, coupling_J, T, ix, iy);
            if (exp(-action_change) > metropolis)
            {
                // accept
                spin[ix][iy] = next_spin; // flip
            }
            else
            {
                // reject
            }
            if (iter > 1000 * nx * ny && (iter + 1) % nskip == 0)
            {
                double m2 = calc_squared_magnetization(spin);
                total_m2 += m2;
                total_m4 += m2*m2;
                count++;
            }
        }
        std::cout << count << std::endl;
        total_m2 /= count;
        total_m4 /= count;

        double binder_ratio = total_m4 / total_m2 / total_m2;
        std::cout << std::fixed << std::setprecision(4)
                  << T << "   "
                  << binder_ratio << "   "
                  << std::endl;
        outputfile << std::fixed << std::setprecision(4)
                   << T << "   "
                   << binder_ratio << "   "
                   << std::endl;
    }
    outputfile.close();
    return 0;
}