#include "gray.h"
#include <limits>
#include <algorithm> // for std::min/std::max

struct Individual {
    std::vector<unsigned int> geno_int; // Gray-coded chromosome (N integers)
    std::vector<double> x;              // Decoded real-valued solution (N doubles)
    double fitness;                     // Function value (lower is better)
};

int f_calls_gray = 0;

RunResult run_single_gray(int bits_per_dim, std::mt19937 &rng, int func_id) 
{
    const double domain_min = (func_id == 1) ? DOMAIN_MIN_FUNC_ONE : DOMAIN_MIN_FUNC_TWO;
    const double domain_max = (func_id == 1) ? DOMAIN_MAX_FUNC_ONE : DOMAIN_MAX_FUNC_TWO;
    const unsigned int maxInt = (bits_per_dim == 32) ? 0xFFFFFFFFu : ((1u << bits_per_dim) - 1u);
    
    // Some distributions that were unused in the original were removed for clarity.
    std::uniform_real_distribution<> dis_prob(0.0, 1.0);
    std::uniform_int_distribution<> dis_pop(0, POPULATION_SIZE - 1);

    std::vector<double> running_best(MAX_F_CALLS, std::numeric_limits<double>::max());
    double overall_best_F = std::numeric_limits<double>::max();
    f_calls_gray = 0;

    // --- 1. INITIALIZATION ---
    std::vector<Individual> population(POPULATION_SIZE);
    std::vector<Individual> next_population;
    next_population.reserve(POPULATION_SIZE);

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population[i].geno_int.resize(N);
        population[i].x.resize(N);
        double sumsq = 0.0, sumcos = 0.0;
        
        for (int j = 0; j < N; ++j) {
            // Randomly initialize the integer value, then convert to Gray code
            std::uniform_int_distribution<unsigned int> dis_val(0, maxInt);
            unsigned int binary_val = dis_val(rng);
            population[i].geno_int[j] = gray_from_binary(binary_val);
            
            unsigned int decoded_int = binary_from_gray(population[i].geno_int[j]);
            population[i].x[j] = int_to_real(decoded_int, bits_per_dim, domain_min, domain_max);
        }
        
        population[i].fitness = calculate_fitness(population[i].x, func_id, N, sumsq, sumcos);
        ++f_calls_gray;

        // Update overall best immediately
        if (population[i].fitness < overall_best_F) {
            overall_best_F = population[i].fitness;
        }

        // Track best-so-far after each function call
        if (f_calls_gray - 1 < MAX_F_CALLS) {
            running_best[f_calls_gray - 1] = overall_best_F;
        }
    }
    
    // --- 2. MAIN GENERATIONAL LOOP ---
    for (int gen = 0; gen < MAX_GENERATIONS && f_calls_gray < MAX_F_CALLS; ++gen) {
        
        // --- Selection (Tournament Selection - size 5) ---
        auto select_parent = [&]() -> const Individual& {
            const Individual* best_ind = nullptr;
            for (int i = 0; i < 5; ++i) { // Tournament size 5
                int idx = dis_pop(rng);
                if (best_ind == nullptr || population[idx].fitness < best_ind->fitness) {
                    best_ind = &population[idx];
                }
            }
            return *best_ind;
        };

        // --- Elitism (Carry the best individual to the next generation) ---
        next_population.clear();
        next_population.reserve(POPULATION_SIZE);
        const Individual* elite = &population[0];
        for (const auto& ind : population) {
            if (ind.fitness < elite->fitness) {
                elite = &ind;
            }
        }
        next_population.push_back(*elite);

        // --- Crossover & Mutation ---
        while (next_population.size() < (size_t)POPULATION_SIZE && f_calls_gray < MAX_F_CALLS) {
            
            // Select two parents
            const Individual& parent1 = select_parent();
            const Individual& parent2 = select_parent();
            
            Individual offspring1 = parent1;
            Individual offspring2 = parent2;

            // Determine whether to do crossover
            if (dis_prob(rng) < CROSSOVER_RATE) {
                int total_bits = N * bits_per_dim;
                // If there is only 1 or fewer total bits, skip crossover to avoid invalid distribution
                if (total_bits > 1) {
                    std::uniform_int_distribution<> dis_crossover_point(1, total_bits - 1);
                    int crossover_point = dis_crossover_point(rng);

                    // Single-point crossover across flattened genome.
                    // Note: bit index order is d*bits_per_dim + b (b = LSB..MSB for each dim).
                    for (int d = 0; d < N; ++d) {
                        for (int b = 0; b < bits_per_dim; ++b) {
                            int current_bit = d * bits_per_dim + b;
                            if (current_bit >= crossover_point) {
                                unsigned int mask = 1u << b;

                                bool bit1 = (parent1.geno_int[d] & mask) != 0;
                                bool bit2 = (parent2.geno_int[d] & mask) != 0;

                                if (bit2) offspring1.geno_int[d] |= mask;
                                else offspring1.geno_int[d] &= ~mask;

                                if (bit1) offspring2.geno_int[d] |= mask;
                                else offspring2.geno_int[d] &= ~mask;
                            }
                        }
                    }
                }
                // else: no valid crossover point -> offspring remain copies of parents
            }
            
            // --- Mutation (Bit-Flip) ---
            auto mutate = [&](Individual& ind) {
                for (int d = 0; d < N; ++d) {
                    for (int b = 0; b < bits_per_dim; ++b) {
                        if (dis_prob(rng) < MUTATION_RATE) {
                            ind.geno_int[d] ^= (1u << b);
                        }
                    }
                }
            };
            
            mutate(offspring1);
            mutate(offspring2);

            // --- Evaluation and Replacement ---
            auto evaluate_and_add = [&](Individual& ind) -> bool {
                if (f_calls_gray >= MAX_F_CALLS) return false; // defensive
                
                double sumsq = 0.0, sumcos = 0.0;
                for (int d = 0; d < N; ++d) {
                    unsigned int decoded_int = binary_from_gray(ind.geno_int[d]);
                    ind.x[d] = int_to_real(decoded_int, bits_per_dim, domain_min, domain_max);
                }
                ind.fitness = calculate_fitness(ind.x, func_id, N, sumsq, sumcos);
                ++f_calls_gray;

                if (ind.fitness < overall_best_F) {
                    overall_best_F = ind.fitness;
                }
                
                // Track best-so-far after each function call
                if (f_calls_gray - 1 < MAX_F_CALLS) {
                    running_best[f_calls_gray - 1] = overall_best_F;
                } else {
                    // if we've exceeded MAX_F_CALLS, do not add to population
                    return false;
                }
                next_population.push_back(std::move(ind));
                return true;
            };

            // Add offspring to the next generation while respecting population size and MAX_F_CALLS
            if (next_population.size() < (size_t)POPULATION_SIZE) {
                if (!evaluate_and_add(offspring1)) break;
            }
            if (next_population.size() < (size_t)POPULATION_SIZE) {
                if (!evaluate_and_add(offspring2)) break;
            }
        }

        // Move next_population into population for next generation
        population = std::move(next_population);
        // Note: next_generation is already reserved and will be reused in next loop iteration.
    }
    
    // Pad remaining function calls with the final best value if the GA finished early
    for (int i = f_calls_gray; i < MAX_F_CALLS; ++i) {
        running_best[i] = overall_best_F;
    }

    return RunResult{running_best, overall_best_F};
}
