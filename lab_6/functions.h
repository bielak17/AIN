#pragma once

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>


const int BITS_PER_DIM = 16;                    // bits per dimension
const double DOMAIN_MIN_FUNC_ONE = -30;        // domain for function 1
const double DOMAIN_MAX_FUNC_ONE = 30;
const double DOMAIN_MIN_FUNC_TWO = -100;     // domain for function 2
const double DOMAIN_MAX_FUNC_TWO = 100;
const double DOMAIN_MIN_FUNC_THREE = -10.24;        // domain for function 1
const double DOMAIN_MAX_FUNC_THREE = 10.24;
const int MAX_F_CALLS = 10000;                  // evaluations per run
const int REPEATS = 100;                        // number of runs
const int N[] = {5,15,30};                      // dimension

const int POPULATION_SIZE = 50;
const int MAX_GENERATIONS = 200;
const double CROSSOVER_RATE = 0.9;
const double MUTATION_RATE = 0.03;

unsigned int binary_from_gray(unsigned int);

unsigned int gray_from_binary(unsigned int);

double int_to_real(unsigned int , int , double , double );

unsigned int real_to_int(double , int , double , double );

double function_one_from_sum_sq(double);

double function_two_from_sums(double, double, int);

double calculate_fitness(const std::vector<double>& , int , int , double& , double& );

double f_rosenbrock(const std::vector<double>& x);

double f_salomon(const std::vector<double>& x);

double f_whitley(const std::vector<double>& x);

void save_2d_function_grid(int func_id, const std::string& filename);

struct RunResult 
{
    std::vector<double> running_best;
    double final_best;
};
