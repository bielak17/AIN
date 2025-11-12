#include "functions.h"



unsigned int binary_from_gray(unsigned int g) 
{
    unsigned int b = 0;
    for (; g; g >>= 1) b ^= g;
    return b;
}

unsigned int gray_from_binary(unsigned int b) 
{
    return b ^ (b >> 1);
}

double int_to_real(unsigned int val, int bits, double domain_min, double domain_max) 
{
    unsigned int maxInt = (bits == 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
    double t = static_cast<double>(val) / static_cast<double>(maxInt);
    return domain_min + t * (domain_max - domain_min);
}

unsigned int real_to_int(double x, int bits, double domain_min, double domain_max)
{
    unsigned int maxInt = (bits == 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
    double t = (x - domain_min) / (domain_max - domain_min);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    return static_cast<unsigned int>(std::round(t * maxInt));
}

double function_one_from_sum_sq(double sumsq) 
{
    double denom = 1.0 + sumsq;
    double inner = -5.0 / denom;                       // -5 / (1 + sum x_i^2)
    double exp_val = std::exp(inner);                  // e^{-5/(1+sum x_i^2)}
    double sin_e = std::sin(exp_val);
    double cos_e = std::cos(exp_val);   
    double cot_e = cos_e / sin_e;
    return inner + std::sin(cot_e);                    // -5/(1+sum) + sin(cot(exp_val))
}

double function_two_from_sums(double sum_sq, double sum_cos, int n) 
{
    double root = std::sqrt(sum_sq / n);                                                            // sqrt(sum x_i^2 / n)
    double value = -20.0 * std::exp(-0.2 * root) - std::exp(sum_cos / n) + 20.0 + std::exp(1.0);    // Ackley's function
    return value;
}

double calculate_fitness(const std::vector<double>& x, int func_id, int n, double& sumsq, double& sumcos)
{
    sumsq = 0.0;
    sumcos = 0.0;

    for (int i = 0; i < n; ++i) {
        double xi = x[i];
        sumsq += xi * xi;
        sumcos += std::cos(2.0 * M_PI * xi);
    }

    double fitness = 0.0;

    if (func_id == 1) {
        // Custom benchmark function one
        fitness = function_one_from_sum_sq(sumsq);
    } 
    else if (func_id == 2) {
        // Ackley-like function two
        fitness = function_two_from_sums(sumsq, sumcos, n);
    } 
    else {
        fitness = 1e9;
    }

    return fitness;
}
