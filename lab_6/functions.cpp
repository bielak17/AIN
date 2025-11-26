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

double f_rosenbrock(const std::vector<double>& x) {
    // Generalised Rosenbrock: sum_{i=0..n-2} [100*(x[i+1]-x[i]^2)^2 + (x[i]-1)^2]
    double sum = 0.0;
    if (x.size() < 2) return 0.0;
    for (size_t i = 0; i + 1 < x.size(); ++i) {
        double t1 = x[i+1] - x[i] * x[i];
        double t2 = x[i] - 1.0;
        sum += 100.0 * t1 * t1 + t2 * t2;
    }
    return sum;
}

double f_salomon(const std::vector<double>& x) {
    // Salomon: 1 - cos(2*pi*sqrt(sum(x_i^2))) + 0.1 * sqrt(sum(x_i^2))
    double sumsq = 0.0;
    for (double xi : x) sumsq += xi * xi;
    double r = std::sqrt(sumsq);
    return 1.0 - std::cos(2.0 * M_PI * r) + 0.1 * r;
}

double f_whitley(const std::vector<double>& x) {
    // Whitley function:
    // f(x) = sum_{i=1..n} sum_{j=1..n} [ ((100*(x_i^2 - x_j)^2 + (1 - x_j)^2)^2)/4000 - cos(100*(x_i^2 - x_j)^2 + (1 - x_j)^2) + 1 ]
    double sum = 0.0;
    int n = static_cast<int>(x.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double a = x[i] - x[j]*x[j];
            double b = 1.0 - x[j];
            double inner = 100.0 * (a * a) + (b * b);
            double term = (inner * inner) / 4000.0 - std::cos(inner) + 1.0;
            sum += term;
        }
    }
    return sum;
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
        fitness = f_rosenbrock(x);
    } 
    else if (func_id == 2) {
        fitness = f_salomon(x);
    } 
    else if (func_id == 3) {
        fitness = f_whitley(x);
    }
    else {
        fitness = 1e9;
    }

    return fitness;
}


void save_2d_function_grid(int func_id, const std::string& filename)
{
    double minX, maxX;
    double minY, maxY;
    switch (func_id)
    {
        case 1: 
            minX = minY = -5.0;
            maxX = maxY = 5.0;
            break;

        case 2:  
            minX = minY = -10.0;
            maxX = maxY = 10.0;
            break;

        case 3:  
            minX = minY = -5.0;
            maxX = maxY = 5.0;
            break;

        default:
            std::cerr << "Invalid func_id in save_2d_function_grid\n";
            return;
    }

    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Could not open file " << filename << "\n";
        return;
    }

    file << std::fixed << std::setprecision(6);
    file << "x,y,z\n";

    std::vector<double> point(2);

    for (double x = minX; x <= maxX; x += 0.1)
    {
        for (double y = minY; y <= maxY; y += 0.1)
        {
            point[0] = x;
            point[1] = y;

            double dummy1, dummy2;
            double z = calculate_fitness(point, func_id, 2, dummy1, dummy2);

            file << x << "," << y << "," << z << "\n";
        }
    }

    file.close();
    std::cout << "Saved 2D grid to " << filename << "\n";
}