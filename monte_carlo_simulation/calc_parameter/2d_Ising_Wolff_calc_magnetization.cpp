#include <iostream>
#include <cmath>
#include <fstream>
const long int niter = 1000000;
const int L = 128;
const int nx = L; // number of sites along x-direction
const int ny = L; // number of sites along y-direction
const double coupling_J = 1.0;
const int nconf = 80;
const double t_start = 1.8;
const int nskip = 10;

double calc_action(const int spin[nx][ny], const double coupling_J, const double temperature)
{
    double action = 0e0;
    int sum = 0;
    for (int ix = 0; ix != nx; ix++)
    {
        int ixp1 = (ix + 1) % nx;
        for (int iy = 0; iy != ny; iy++)
        {
            int iyp1 = (iy + 1) % ny;
            sum = sum + spin[ix][iy] * spin[ixp1][iy] + spin[ix][iy] * spin[ix][iyp1];
        }
    }
    action = sum * coupling_J / temperature * (-1e0);

    return action;
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
    if(total_spin < 0)
    {
        total_spin = -total_spin;
    }
    return total_spin;
}

double calc_magnetization(const int total_spin, const int total_sampling)
{
    double magnetization = 0;
    magnetization = total_spin / total_sampling;
    if (magnetization < 0)
    {
        magnetization = -magnetization;
    }
    return magnetization;
}

int make_cluster(const int spin[nx][ny], const double coupling_J,
                 const double temperature, int &n_cluster, int (&i_cluster)[nx * ny][2])
{
    int in_or_out[nx][ny];
    for (int ix = 0; ix != nx; ix++)
    {
        for (int iy = 0; iy != ny; iy++)
        {
            in_or_out[ix][iy] = 1;
        }
    }
    // in_or_out[ix][iy] = 1 -> not in the cluster; 0 -> in the cluster.
    double rand_site = (double)rand() / RAND_MAX;
    rand_site = rand_site * nx * ny;
    int ix = (int)rand_site / ny;
    int iy = (int)rand_site % ny;
    in_or_out[ix][iy] = 0;
    i_cluster[0][0] = ix;
    i_cluster[0][1] = iy;
    int spin_cluster = spin[ix][iy];
    n_cluster = 1;
    double probability = 1e0 - exp(-2e0 * coupling_J / temperature);
    int k = 0;
    while (k < n_cluster)
    {
        ix = i_cluster[k][0];
        iy = i_cluster[k][1];
        int ixp1 = (ix + 1) % nx;
        int iyp1 = (iy + 1) % ny;
        int ixm1 = (ix - 1 + nx) % nx;
        int iym1 = (iy - 1 + ny) % ny;

        if (spin[ixp1][iy] == spin_cluster)
        {
            if (in_or_out[ixp1][iy] == 1)
            {
                if ((double)rand() / RAND_MAX < probability)
                {
                    i_cluster[n_cluster][0] = ixp1;
                    i_cluster[n_cluster][1] = iy;
                    n_cluster = n_cluster + 1;
                    in_or_out[ixp1][iy] = 0;
                }
            }
        }
        if (spin[ix][iyp1] == spin_cluster)
        {
            if (in_or_out[ix][iyp1] == 1)
            {
                if ((double)rand() / RAND_MAX < probability)
                {
                    i_cluster[n_cluster][0] = ix;
                    i_cluster[n_cluster][1] = iyp1;
                    n_cluster = n_cluster + 1;
                    in_or_out[ix][iyp1] = 0;
                }
            }
        }
        if (spin[ixm1][iy] == spin_cluster)
        {
            if (in_or_out[ixm1][iy] == 1)
            {
                if ((double)rand() / RAND_MAX < probability)
                {
                    i_cluster[n_cluster][0] = ixm1;
                    i_cluster[n_cluster][1] = iy;
                    n_cluster = n_cluster + 1;
                    in_or_out[ixm1][iy] = 0;
                }
            }
        }
        if (spin[ix][iym1] == spin_cluster)
        {
            if (in_or_out[ix][iym1] == 1)
            {
                if ((double)rand() / RAND_MAX < probability)
                {
                    i_cluster[n_cluster][0] = ix;
                    i_cluster[n_cluster][1] = iym1;
                    n_cluster = n_cluster + 1;
                    in_or_out[ix][iym1] = 0;
                }
            }
        }
        k = k + 1;
    }
    return spin_cluster;
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
    std::ofstream outputfile("output/2d_Ising_L" + std::to_string(L) + "_magnetization_wolff.txt");
    for (int conf = 0; conf < nconf + 1; conf++)
    {
        double T = temperature[conf];
        int count = 0;
        int total_spin_sum = 0;
        int spin[nx][ny];
        srand((unsigned)time(NULL));
        for (int ix = 0; ix != nx; ix++)
        {
            for (int iy = 0; iy != ny; iy++)
            {
                spin[ix][iy] = 1;
            }
        }

        for (long int iter = 0; iter != niter; iter++)
        {
            int n_cluster;
            int i_cluster[nx * ny][2];
            int spin_cluster = make_cluster(spin, coupling_J, T, n_cluster, i_cluster);
            double metropolis = (double)rand() / RAND_MAX;
            for (int k = 0; k != n_cluster; k++)
            {
                int ix = i_cluster[k][0];
                int iy = i_cluster[k][1];
                spin[ix][iy] = -spin[ix][iy];
            }

            if (iter > 100000 && (iter + 1) % nskip == 0)
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
                  << m << std::endl;
        outputfile << std::fixed << std::setprecision(4)
                   << T << "   "
                   << m << std::endl;
    }
    outputfile.close();
    return 0;
}
