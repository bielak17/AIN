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

const int BITS_PER_DIM = 16;
const double DOMAIN_MIN = -10.0;
const double DOMAIN_MAX = 10.0;
const int MAX_F_CALLS = 10000;      // evaluations per run
const int REPEATS = 100;            // number of runs
const int NS[] = {2, 5, 10};        // array of dimensions

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


inline double int_to_real(unsigned int val, int bits) 
{
    unsigned int maxInt = (bits == 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
    double t = static_cast<double>(val) / static_cast<double>(maxInt);
    return DOMAIN_MIN + t * (DOMAIN_MAX - DOMAIN_MIN);
}

inline unsigned int real_to_int(double x, int bits)
{
    unsigned int maxInt = (bits == 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
    double t = (x - DOMAIN_MIN) / (DOMAIN_MAX - DOMAIN_MIN);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    return static_cast<unsigned int>(std::round(t * maxInt));
}


struct RunResult 
{
    std::vector<double> running_best;
    double final_best;
};


RunResult run_single(int n, int bits_per_dim, std::mt19937 &rng) 
{
    int total_bits = n * bits_per_dim;

    std::vector<unsigned int> geno_int(n, 0);
    std::vector<unsigned int> decoded_int(n, 0);
    std::vector<double> x(n, 0.0);
    std::vector<double> x2(n, 0.0);

    unsigned int maxInt = (bits_per_dim == 32) ? 0xFFFFFFFFu : ((1u << bits_per_dim) - 1u);

    for (int i = 0; i < n; ++i) 
    {
        double start = DOMAIN_MAX;
        unsigned int intval = real_to_int(start, bits_per_dim);
        geno_int[i] = gray_from_binary(intval);
        decoded_int[i] = intval;
        x[i] = int_to_real(decoded_int[i], bits_per_dim);
        x2[i] = x[i] * x[i];
    }

    double F = 0.0;
    for (int i = 0; i < n; ++i) F += x2[i];
    int f_calls = 1;

    std::vector<double> running(MAX_F_CALLS, F);

    std::vector<std::pair<int, int>> neighbor_list;
    neighbor_list.reserve(total_bits);
    for (int d = 0; d < n; ++d)
        for (int b = 0; b < bits_per_dim; ++b)
            neighbor_list.emplace_back(d, b);

    bool improved = true;
    while (f_calls < MAX_F_CALLS && improved) 
    {
        improved = false;
        std::shuffle(neighbor_list.begin(), neighbor_list.end(), rng);

        for (auto &nb : neighbor_list) {
            if (f_calls >= MAX_F_CALLS) break;

            int dim = nb.first;
            int bitpos = nb.second;

            unsigned int gcur = geno_int[dim];
            unsigned int gnew = gcur ^ (1u << bitpos);

            unsigned int bnew = binary_from_gray(gnew);
            bnew &= maxInt;

            double x_old = x[dim];
            double x_new = int_to_real(bnew, bits_per_dim);
            double F_new = F - x_old * x_old + x_new * x_new;

            if (f_calls < MAX_F_CALLS) 
            {
                double prev_best = running[f_calls - 1];
                running[f_calls] = std::min(prev_best, F_new);
            }
            f_calls++;

            if (F_new < F - 1e-15) 
            {
                geno_int[dim] = gnew;
                decoded_int[dim] = bnew;
                x[dim] = x_new;
                x2[dim] = x_new * x_new;
                F = F_new;
                improved = true;
                if (f_calls - 1 < MAX_F_CALLS) 
                {
                    running[f_calls - 1] = std::min(running[f_calls - 1], F);
                }
                break;
            }
        }
    }

    for (int k = 1; k < MAX_F_CALLS; ++k)
        if (running[k] > running[k - 1]) running[k] = running[k - 1];

    RunResult rr;
    rr.running_best = std::move(running);
    rr.final_best = F;
    return rr;
}

int main() 
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::random_device rd;

    for (int n : NS) 
    {
        std::vector<std::vector<double>> all_runs(REPEATS, std::vector<double>(MAX_F_CALLS));

        for (int run = 0; run < REPEATS; ++run) 
        {
            std::mt19937 run_rng(rd() + run * 7919 + n * 101);
            RunResult rr = run_single(n, BITS_PER_DIM, run_rng);
            all_runs[run] = std::move(rr.running_best);
        }

        std::ostringstream fname;
        fname << "runs_n" << n << ".csv";
        std::ofstream of(fname.str());

        for (int r = 0; r < REPEATS; ++r) 
        {
            of << "run_" << (r + 1);
            if (r < REPEATS - 1) of << ",";
        }
        of << "\n";

        for (int step = 0; step < MAX_F_CALLS; ++step) 
        {
            for (int r = 0; r < REPEATS; ++r) 
            {
                of << all_runs[r][step];
                if (r < REPEATS - 1) of << ",";
            }
            of << "\n";
        }

        of.close();
    }

    std::cout << "All experiments completed. CSV files saved in current folder.\n";
    system("pause");
}
