#include <iostream>
#include <random>
#include <math.h>
// #include <pybind11/pybind11.h>

using namespace std;

double get_uniform_distributed_rand(double min_val, double max_val)
{
    /***** 一様乱数生成 *****/
    random_device seed;
    mt19937 engine(seed()); // メルセンヌ・ツイスター法
    // std::minstd_rand0 engine(seed());    // 線形合同法
    // std::ranlux24_base engine(seed());   // キャリー付き減算法

    uniform_real_distribution<double> get_rand(min_val, max_val);

    return get_rand(engine);
}

double get_normal_distributed_rand(double mu, double sigma)
{
    /***** 一次元正規乱数生成 *****/
    random_device seed;
    mt19937 engine(seed()); // メルセンヌ・ツイスター法
    // std::minstd_rand0 engine(seed());    // 線形合同法
    // std::ranlux24_base engine(seed());   // キャリー付き減算法

    normal_distribution<double> get_rand(mu, sigma);
    return get_rand(engine);
}
