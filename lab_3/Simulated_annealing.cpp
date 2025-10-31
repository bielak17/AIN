#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>
#include <utility>

const int BITS_PER_DIM = 16;                    // bits per dimension
const double DOMAIN_MIN_FUNC_ONE = -3.0;        // domain for function 1
const double DOMAIN_MAX_FUNC_ONE = 3.0;
const double DOMAIN_MIN_FUNC_TWO = -32.768;     // domain for function 2
const double DOMAIN_MAX_FUNC_TWO = 32.768;
const int MAX_F_CALLS = 10000;                  // evaluations per run
const int REPEATS = 100;                        // number of runs
const int N = 10;                               // dimension

inline unsigned int binary_from_gray(unsigned int g) 
{
    unsigned int b = 0;
    for (; g; g >>= 1) b ^= g;
    return b;
}

inline unsigned int gray_from_binary(unsigned int b) 
{
    return b ^ (b >> 1);
}

inline double int_to_real(unsigned int val, int bits, double domain_min, double domain_max) 
{
    unsigned int maxInt = (bits == 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
    double t = static_cast<double>(val) / static_cast<double>(maxInt);
    return domain_min + t * (domain_max - domain_min);
}

inline unsigned int real_to_int(double x, int bits, double domain_min, double domain_max)
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

struct RunResult 
{
    std::vector<double> running_best;
    double final_best;
};


RunResult run_single_gray(int bits_per_dim, std::mt19937 &rng, int func_id) 
{
    const double domain_min = (func_id == 1) ? DOMAIN_MIN_FUNC_ONE : DOMAIN_MIN_FUNC_TWO;
    const double domain_max = (func_id == 1) ? DOMAIN_MAX_FUNC_ONE : DOMAIN_MAX_FUNC_TWO;
    int total_bits = N * bits_per_dim;
    std::vector<unsigned int> geno_int(N, 0);
    std::vector<unsigned int> decoded_int(N, 0);
    std::vector<double> x(N, 0.0);

    unsigned int maxInt = (bits_per_dim == 32) ? 0xFFFFFFFFu : ((1u << bits_per_dim) - 1u);

    for (int i = 0; i < N; ++i) 
    {
        double start = domain_max;
        unsigned int intval = real_to_int(start, bits_per_dim, domain_min, domain_max);
        geno_int[i] = gray_from_binary(intval);
        decoded_int[i] = intval;
        x[i] = int_to_real(decoded_int[i], bits_per_dim, domain_min, domain_max);
    }
    double sumsq = 0.0;
    for (int i=0; i<x.size(); i++) sumsq += x[i] * x[i];
    double sumcos = 0.0;
    for (int i=0; i<x.size(); i++) sumcos += std::cos(2.0 * M_PI * x[i]);
    double F = (func_id == 1) ? function_one_from_sum_sq(sumsq) : function_two_from_sums(sumsq, sumcos, N);
    int f_calls = 1;
    double best_F = F;

    std::vector<double> running(MAX_F_CALLS, F);

    //hardcoded initial temperature for testing
    double T1 = 0.1;
    double TN = 0.0001;
    int L = 100;
    int n = 100;
    int k=0;
    double Tk = T1;

    std::vector<std::pair<int, int>> neighbor_list;
    neighbor_list.reserve(total_bits);
    for (int d = 0; d < N; ++d)
        for (int b = 0; b < bits_per_dim; ++b)
            neighbor_list.emplace_back(d, b);

    while (f_calls < MAX_F_CALLS && Tk > TN) 
    {
        std::shuffle(neighbor_list.begin(), neighbor_list.end(), rng);
        for (auto &nb : neighbor_list) {
            if (f_calls >= MAX_F_CALLS) break;
            
            // creating neighbor by flipping one bit
            int dim = nb.first;
            int bitpos = nb.second;

            unsigned int gcur = geno_int[dim];
            unsigned int gnew = gcur ^ (1u << bitpos);

            unsigned int bnew = binary_from_gray(gnew);
            bnew &= maxInt;

            double x_old = x[dim];
            double x_new = int_to_real(bnew, bits_per_dim, domain_min, domain_max);
            double new_sum_sq = sumsq - x_old * x_old + x_new * x_new;
            double new_sum_cos = sumcos - std::cos(2.0 * M_PI * x_old) + std::cos(2.0 * M_PI * x_new);
            double F_new = (func_id == 1) ? function_one_from_sum_sq(new_sum_sq) : function_two_from_sums(new_sum_sq, new_sum_cos, N);

            //saving current best in running
            if (f_calls < MAX_F_CALLS) 
            {
                double prev_best = running[f_calls - 1];
                running[f_calls] = std::min(prev_best, F_new);
            }
            f_calls++;

            // accepting improvement and updating current best solution
            if (F_new < F - 1e-15) 
            {
                geno_int[dim] = gnew;
                decoded_int[dim] = bnew;
                x[dim] = x_new;
                F = F_new;
                best_F = std::min(best_F, F);
                sumsq = new_sum_sq;
                sumcos = new_sum_cos;
                if (f_calls - 1 < MAX_F_CALLS) 
                {
                    running[f_calls - 1] = std::min(running[f_calls - 1], F);
                }
            }
            // else if for probabilistic acceptance here
            else if (std::exp((F - F_new) / Tk) > std::uniform_real_distribution<>(0.0, 1.0)(rng)) 
            {
                geno_int[dim] = gnew;
                decoded_int[dim] = bnew;
                x[dim] = x_new;
                F = F_new;
                best_F = std::min(best_F, F);
                sumsq = new_sum_sq;
                sumcos = new_sum_cos;
                if (f_calls - 1 < MAX_F_CALLS) 
                {
                    running[f_calls - 1] = std::min(running[f_calls - 1], F);
                }
            }
        }
        k++;
        Tk = (T1-TN) / (1+std::exp(0.3*(k-n/2)))+TN;      // cooling
        // debug print
        //std::cout << "Func " << func_id << " Temperature: " << T << ", Best F: " << best_F << std::endl;
    }

    // Change to ensure running best is non-increasing
    for (int k = 1; k < MAX_F_CALLS; ++k)
        if (running[k] > running[k - 1]) running[k] = running[k - 1];

    // creating and returning best result
    RunResult rr;
    rr.running_best = std::move(running);
    rr.final_best = best_F;
    return rr;
}

RunResult run_single_real(int bits_per_dim, std::mt19937 &rng, int func_id) 
{
    const double domain_min = (func_id == 1) ? DOMAIN_MIN_FUNC_ONE : DOMAIN_MIN_FUNC_TWO;
    const double domain_max = (func_id == 1) ? DOMAIN_MAX_FUNC_ONE : DOMAIN_MAX_FUNC_TWO;
    int total_bits = N * bits_per_dim;

    std::vector<double> x(N, 0.0);

    unsigned int maxInt = (bits_per_dim == 32) ? 0xFFFFFFFFu : ((1u << bits_per_dim) - 1u);
    for (int i = 0; i < N; ++i) 
    {
        double start = domain_max;
        x[i] = start;
    }
    double sumsq = 0.0;
    for (int i=0; i<x.size(); i++) sumsq += x[i] * x[i];
    double sumcos = 0.0;
    for (int i=0; i<x.size(); i++) sumcos += std::cos(2.0 * M_PI * x[i]);
    double F = (func_id == 1) ? function_one_from_sum_sq(sumsq) : function_two_from_sums(sumsq, sumcos, N);
    int f_calls = 1;
    double best_F = F;

    std::vector<double> running(MAX_F_CALLS, F);
    
    //initial temperature
    double T1 = 0.1;
    double TN = 0.0001;
    int L = 100;
    int n = 100;
    int k=0;
    double Tk = T1;

    while (f_calls < MAX_F_CALLS && Tk > TN) 
    {
        for (int i = 0; i < L; ++i) {
            if (f_calls >= MAX_F_CALLS) break;

            // creating neighbor by adding random number from this generator: N(0,1)
            int dim = rng() % N;
            std::normal_distribution<double> nd(0.0, 1.0);
            double x_new = x[dim] + nd(rng);
            double x_old = x[dim];
            if (x_new < domain_min) x_new = domain_min;
            if (x_new > domain_max) x_new = domain_max;
            double new_sum_sq = sumsq - x_old * x_old + x_new * x_new;
            double new_sum_cos = sumcos - std::cos(2.0 * M_PI * x_old) + std::cos(2.0 * M_PI * x_new);
            double F_new = (func_id == 1) ? function_one_from_sum_sq(new_sum_sq) : function_two_from_sums(new_sum_sq, new_sum_cos, N);

            //saving current best in running
            if (f_calls < MAX_F_CALLS) 
            {
                double prev_best = running[f_calls - 1];
                running[f_calls] = std::min(prev_best, F_new);
            }
            f_calls++;

            // accepting improvement and updating current best solution
            if (F_new < F - 1e-15) 
            {
                x[dim] = x_new;
                F = F_new;
                best_F = std::min(best_F, F);
                sumsq = new_sum_sq;
                sumcos = new_sum_cos;
                if (f_calls - 1 < MAX_F_CALLS) 
                {
                    running[f_calls - 1] = std::min(running[f_calls - 1], F);
                }
            }
            // else if for probabilistic acceptance here
            else if (std::exp((F - F_new) / Tk) > std::uniform_real_distribution<>(0.0, 1.0)(rng)) 
            {
                x[dim] = x_new;
                F = F_new;
                best_F = std::min(best_F, F);
                sumsq = new_sum_sq;
                sumcos = new_sum_cos;
                if (f_calls - 1 < MAX_F_CALLS) 
                {
                    running[f_calls - 1] = std::min(running[f_calls - 1], F);
                }
            }
        }
        k++;
        Tk = (T1-TN) / (1+std::exp(0.3*(k-n/2)))+TN;      // cooling  function e) from lecture
        // debug print
        //std::cout << "Func " << func_id << " Temperature: " << T << ", Best F: " << best_F << std::endl;
    }

    // Change to ensure running best is non-increasing
    for (int k = 1; k < MAX_F_CALLS; ++k)
        if (running[k] > running[k - 1]) running[k] = running[k - 1];

    // creating and returning best result
    RunResult rr;
    rr.running_best = std::move(running);
    rr.final_best = best_F;
    return rr;
}

int main() 
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Running experiments for Gray encoding
    std::random_device rd;
    std::vector<std::vector<double>> all_runs(REPEATS, std::vector<double>(MAX_F_CALLS));
    for (int func_id = 1; func_id <= 2; ++func_id)
    {
        for (int run = 0; run < REPEATS; ++run) 
        {
            std::mt19937 run_rng(rd() + run * 7919 + N * 101);
            RunResult rr = run_single_gray(BITS_PER_DIM, run_rng, func_id);
            all_runs[run] = std::move(rr.running_best);
        }

        // Saving results to CSV
        std::ostringstream fname;
        fname << "results_gray_func_" << func_id << ".csv";
        std::ofstream of(fname.str());

        for (int r = 0; r < REPEATS; ++r) 
        {
            of << "run_" << (r + 1);
            if (r < REPEATS - 1) of << ";";
        }
        of << "\n";

        for (int step = 0; step < MAX_F_CALLS; ++step) 
        {
            for (int r = 0; r < REPEATS; ++r) 
            {
                of << all_runs[r][step];
                if (r < REPEATS - 1) of << ";";
            }
            of << "\n";
        }
        of.close();
    }

    // Running experiments for Real encoding
    for (int func_id = 1; func_id <= 2; ++func_id)
    {
        for (int run = 0; run < REPEATS; ++run) 
        {
            std::mt19937 run_rng(rd() + run * 7919 + N * 101);
            RunResult rr = run_single_real(BITS_PER_DIM, run_rng, func_id);
            all_runs[run] = std::move(rr.running_best);
        }

        // Saving results to CSV
        std::ostringstream fname;
        fname << "results_real_func_" << func_id << ".csv";
        std::ofstream of(fname.str());

        for (int r = 0; r < REPEATS; ++r) 
        {
            of << "run_" << (r + 1);
            if (r < REPEATS - 1) of << ";";
        }
        of << "\n";

        for (int step = 0; step < MAX_F_CALLS; ++step) 
        {
            for (int r = 0; r < REPEATS; ++r) 
            {
                of << all_runs[r][step];
                if (r < REPEATS - 1) of << ";";
            }
            of << "\n";
        }
        of.close();
    }

    std::cout << "All four experiments completed. CSV files saved in current folder.\n";
    return 0;
}