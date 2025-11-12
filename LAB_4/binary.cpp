#include "binary.h"


struct Individual {
    std::vector<unsigned int> geno_int; 
    std::vector<double> x;              // Decoded real-valued solution (N doubles)
    double fitness;                     // Function value (lower is better)
};

// Global variable to track function calls, matching the SA structure
int f_calls_binary = 0;


RunResult run_single_real(int bits_per_dim, std::mt19937 &rng, int func_id) 
{
    const double domain_min = (func_id == 1) ? DOMAIN_MIN_FUNC_ONE : DOMAIN_MIN_FUNC_TWO;
    const double domain_max = (func_id == 1) ? DOMAIN_MAX_FUNC_ONE : DOMAIN_MAX_FUNC_TWO;
    const unsigned int maxInt = (bits_per_dim == 32) ? 0xFFFFFFFFu : ((1u << bits_per_dim) - 1u);
    
    std::uniform_real_distribution<> dis_prob(0.0, 1.0);
    std::uniform_int_distribution<> dis_pop(0, POPULATION_SIZE - 1);
    
    std::vector<double> running_best(MAX_F_CALLS, std::numeric_limits<double>::max());
    double overall_best_F = std::numeric_limits<double>::max();
    f_calls_binary = 0;

    // --- 1. INITIALIZATION ---
    std::vector<Individual> population(POPULATION_SIZE);
    std::vector<Individual> next_population;

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population[i].geno_int.resize(N);
        population[i].x.resize(N);
        double sumsq, sumcos;
        
        for (int j = 0; j < N; ++j) {
            std::uniform_int_distribution<unsigned int> dis_val(0, maxInt);
            unsigned int binary_val = dis_val(rng);
            population[i].geno_int[j] = binary_val;
            
            population[i].x[j] = int_to_real(binary_val, bits_per_dim, domain_min, domain_max);
        }
        
        population[i].fitness = calculate_fitness(population[i].x, func_id, N, sumsq, sumcos);
        f_calls_binary++;
        
        // Update best fitness from the initial population
        if (population[i].fitness < overall_best_F) {
            overall_best_F = population[i].fitness;
        }
    }
    
    if (f_calls_binary - 1 < MAX_F_CALLS) {
        running_best[f_calls_binary - 1] = overall_best_F;
    }
    
    // --- 2. MAIN GENERATIONAL LOOP ---
    for (int gen = 0; gen < MAX_GENERATIONS && f_calls_binary < MAX_F_CALLS; ++gen) {
        
        // --- Selection (Tournament Selection - size 5) ---
        auto select_parent = [&]() -> const Individual& {
            const Individual* best_ind = nullptr;
            for (int i = 0; i < 5; ++i) { 
                int idx = dis_pop(rng);
                if (best_ind == nullptr || population[idx].fitness < best_ind->fitness) {
                    best_ind = &population[idx];
                }
            }
            return *best_ind;
        };

        // --- Elitism (Carry the best individual to the next generation) ---
        next_population.clear();
        const Individual* elite = &population[0];
        for(const auto& ind : population) {
            if (ind.fitness < elite->fitness) {
                elite = &ind;
            }
        }
        next_population.push_back(*elite);

        // --- Crossover & Mutation ---
        while (next_population.size() < POPULATION_SIZE && f_calls_binary < MAX_F_CALLS) {
            
            // Select two parents
            const Individual& parent1 = select_parent();
            const Individual& parent2 = select_parent();
            
            Individual offspring1 = parent1;
            Individual offspring2 = parent2;

            if (dis_prob(rng) < CROSSOVER_RATE) {
                // Single-Point Crossover across the entire N*bits_per_dim genome
                int total_bits = N * bits_per_dim;
                std::uniform_int_distribution<> dis_crossover_point(1, total_bits - 1);
                int crossover_point = dis_crossover_point(rng);

                for (int d = 0; d < N; ++d) {
                    for (int b = 0; b < bits_per_dim; ++b) {
                        int current_bit = d * bits_per_dim + b;
                        
                        if (current_bit >= crossover_point) {
                            // Swap bits after the crossover point
                            unsigned int mask = 1u << b;
                            
                            // Get bit b from parent1's dim d
                            bool bit1 = (parent1.geno_int[d] & mask) != 0;
                            // Get bit b from parent2's dim d
                            bool bit2 = (parent2.geno_int[d] & mask) != 0;
                            
                            // Set bit b for offspring1
                            if (bit2) offspring1.geno_int[d] |= mask;
                            else offspring1.geno_int[d] &= ~mask;
                            
                            // Set bit b for offspring2
                            if (bit1) offspring2.geno_int[d] |= mask;
                            else offspring2.geno_int[d] &= ~mask;
                        }
                    }
                }
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
            auto evaluate_and_add = [&](Individual& ind) {
                double sumsq, sumcos;
                for (int d = 0; d < N; ++d) {
                    // Decoding is direct from binary
                    ind.x[d] = int_to_real(ind.geno_int[d], bits_per_dim, domain_min, domain_max);
                }
                ind.fitness = calculate_fitness(ind.x, func_id, N, sumsq, sumcos);
                f_calls_binary++;
                
                if (ind.fitness < overall_best_F) {
                    overall_best_F = ind.fitness;
                }
                
                // Track best-so-far
                if (f_calls_binary - 1 < MAX_F_CALLS) {
                    running_best[f_calls_binary - 1] = overall_best_F;
                } else {
                    // Stop if MAX_F_CALLS is exceeded
                    return; 
                }
                next_population.push_back(std::move(ind));
            };

            // Add offspring to the next generation
            if (next_population.size() < POPULATION_SIZE) evaluate_and_add(offspring1);
            if (next_population.size() < POPULATION_SIZE) evaluate_and_add(offspring2);
        }

        population = std::move(next_population);
    }
    
    for (int i = f_calls_binary; i < MAX_F_CALLS; ++i) {
        running_best[i] = overall_best_F;
    }

    return RunResult{running_best, overall_best_F};
}
