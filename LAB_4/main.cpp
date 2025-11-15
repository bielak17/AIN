//compile g++ -fdiagnostics-color=always -g main.cpp gray.cpp real.cpp functions.cpp -o main.exe

#include "gray.h"
#include "real.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <sstream>
#include <string>


int main() 
{
    std::cout << "Starting experiments for Real and Gray encoding...\n"<<std::endl;
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
        std::cout << fname.str() << " file created." << std::endl;
    }

    // Running experiments for Real encoding
    for (int func_id = 1; func_id <= 2; ++func_id)
    {
        for (int run = 0; run < REPEATS; ++run) 
        {
            std::mt19937 run_rng(rd() + run * 7919 + N * 101);
            RunResult rr = run_single_real(run_rng, func_id);
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
        std::cout << fname.str() << " file created." << std::endl;
    }

    std::cout << "All four experiments completed. CSV files saved in current folder.\n";
    return 0;
}
