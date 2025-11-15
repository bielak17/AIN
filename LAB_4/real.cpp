#include "real.h"

struct Individual_real {
    std::vector<double> x;              // Decoded real-valued solution (N doubles)
    double sigma;                       // Mutation step size
    double fitness;                     // Function value (lower is better)
};

// Global variable to track function calls, matching the SA structure
int f_calls_binary = 0;



RunResult run_single_real(std::mt19937 &rng, int func_id) 
{
    const double domain_min = (func_id == 1) ? DOMAIN_MIN_FUNC_ONE : DOMAIN_MIN_FUNC_TWO;
    const double domain_max = (func_id == 1) ? DOMAIN_MAX_FUNC_ONE : DOMAIN_MAX_FUNC_TWO;
        
    std::uniform_real_distribution<> dis_prob(0.0, 1.0);
    std::uniform_int_distribution<> dis_pop(0, POPULATION_SIZE - 1);
    
    std::vector<double> running_best(MAX_F_CALLS, std::numeric_limits<double>::max());
    double overall_best_F = std::numeric_limits<double>::max();
    f_calls_binary = 0;

    // --- 1. INITIALIZATION ---
    std::vector<Individual_real> population(POPULATION_SIZE);
    std::vector<Individual_real> next_population;

    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population[i].x.resize(N);
        double sumsq, sumcos;
        
        for (int j = 0; j < N; ++j) {
            std::uniform_real_distribution<double> init_val(domain_min,domain_max);
            population[i].x[j] = init_val(rng);
        }
        
        population[i].sigma = 0.1 * (domain_max - domain_min);
        population[i].fitness = calculate_fitness(population[i].x, func_id, N, sumsq, sumcos);
        f_calls_binary++;
        
        // Update best fitness from the initial population
        if (population[i].fitness < overall_best_F) {
            overall_best_F = population[i].fitness;
        }

    
        if (f_calls_binary - 1 < MAX_F_CALLS) {
            running_best[f_calls_binary - 1] = overall_best_F;
        }
    }
    
    // --- 2. MAIN GENERATIONAL LOOP ---
    for (int gen = 0; gen < MAX_GENERATIONS && f_calls_binary < MAX_F_CALLS; ++gen) {
        
        // --- Selection (Tournament Selection - size 5) ---
        auto select_parent = [&]() -> const Individual_real& {
            const Individual_real* best_ind = nullptr;
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
        const Individual_real* elite = &population[0];
        for(const auto& ind : population) {
            if (ind.fitness < elite->fitness) {
                elite = &ind;
            }
        }
        next_population.push_back(*elite);

        // --- Crossover & Mutation ---
        while (next_population.size() < POPULATION_SIZE && f_calls_binary < MAX_F_CALLS) {
            
            // Select multiple parents
            const Individual_real& parent1 = select_parent();
            const Individual_real& parent2 = select_parent();
            const Individual_real& parent3 = select_parent();
            
            Individual_real offspring1 = parent1;

            if (dis_prob(rng) < CROSSOVER_RATE) {
                //Classic intermediate real valued recombination - average from all parents
                for (int i = 0; i < N; ++i) {
                    offspring1.x[i] = (parent1.x[i] + parent2.x[i] + parent3.x[i]) / 3.0;
                }
            }
            
            // --- Mutation #1 self-adaptive from lecture ---
            auto mutate = [&](Individual_real& ind) {
                ind.sigma *= std::exp(1/std::sqrt(N) * std::normal_distribution<>(0.0,1.0)(rng));
                if (ind.sigma < 1e-8)
                    ind.sigma = 1e-8;
                for (int i=0; i < N; ++i) {
                    ind.x[i] += ind.sigma * std::normal_distribution<>(0.0,ind.sigma)(rng);
                    if (ind.x[i] < domain_min)
                        ind.x[i] = domain_min;
                    if (ind.x[i] > domain_max)
                        ind.x[i] = domain_max;
                }
            };
            
            mutate(offspring1);

            // --- Evaluation and Replacement ---
            auto evaluate_and_add = [&](Individual_real& ind) {
                double sumsq, sumcos;
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
            if (next_population.size() < POPULATION_SIZE)
                evaluate_and_add(offspring1);
        }
        population = std::move(next_population);
    }
    
    for (int i = f_calls_binary; i < MAX_F_CALLS; ++i) {
        running_best[i] = overall_best_F;
    }

    return RunResult{running_best, overall_best_F};
}
