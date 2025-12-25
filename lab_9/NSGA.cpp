#include "NSGA.h"

void evaluateZDT1(Individual& a)
{
    double f1 = a.values[0];
    double g = 1+9*a.values[1];
    double h = 1-sqrt(f1/g);
    double f2 = g*h;
    a.objectives[0] = f1;
    a.objectives[1] = f2;
}

bool dominates(const vector<double>& a, const vector<double>& b)
{
    bool strictlyBetter = false;
    for (size_t i = 0; i < a.size(); i++)
    {
        if (a[i] > b[i])
            return false;
        if (a[i] < b[i])
            strictlyBetter = true;
    }
    return strictlyBetter;
}

void paretoLayers(vector<Individual>& population)
{
    vector<Individual*> remaining;
    for (auto& p : population)
        remaining.push_back(&p);
    int currentRank = 1;
    while (!remaining.empty())
    {
        //Find current Pareto front
        vector<Individual*> front;
        for (auto* p : remaining)
        {
            bool dominated = false;
            for (auto* q : remaining)
            {
                if (p != q && dominates(q->objectives, p->objectives))
                {
                    dominated = true;
                    break;
                }
            }
            if (!dominated)
                front.push_back(p);
        }
        // Assign rank
        for (auto* p : front)
            p->rank = currentRank;
        // Remove front from remaining
        vector<Individual*> newRemaining;
        for (auto* p : remaining)
        {
            if (find(front.begin(), front.end(), p) == front.end())
                newRemaining.push_back(p);
        }
        // Update remaining population and apply next rank
        remaining = newRemaining;
        currentRank++;
    }
}

void calculateCrowdingDistances(vector<Individual>& population)
{
    // Set all dist to 0
    for (auto& ind : population)
    {
        ind.crowding_distance = 0.0;
    }
    // Find max rank
    int maxRank = 1;
    for (const auto& ind : population)
    {
        if (ind.rank > maxRank)
            maxRank = ind.rank;
    }
    for (int rank = 1; rank <= maxRank; rank++)
    {
        // Collect all individuals from this front
        vector<Individual*> front;
        for (auto& ind : population)
        {
            if (ind.rank == rank)
                front.push_back(&ind);
        }
        // Calculate crowding distance for each objective
        for (int m = 0; m < 2; m++)
        {
            // Sort front by objective m (1 or 2)
            sort(front.begin(), front.end(),
                [m](const Individual* a, const Individual* b)
                {
                    return a->objectives[m] < b->objectives[m];
                });
            // set infinite distance to boundary points
            front[0]->crowding_distance = numeric_limits<double>::infinity();
            front.back()->crowding_distance = numeric_limits<double>::infinity();
            double f_min = front[0]->objectives[m];
            double f_max = front.back()->objectives[m];
            // Calculate crowding distance for other points
            for (size_t i = 1; i < front.size() - 1; i++)
            {
                if (f_max - f_min == 0) 
                    continue; // Avoid division by zero
                front[i]->crowding_distance += (front[i + 1]->objectives[m] - front[i - 1]->objectives[m]) / (f_max - f_min);
            }
        }
    }
}

vector<Individual> NSGA_II(mt19937 &rng, int dims)
{
    int f_calls = 0;
    // --- 1. INITIALIZATION ---
    uniform_real_distribution<double> init_val(DOMAIN_MIN,DOMAIN_MAX);
    uniform_int_distribution<int> pop_dist(0, POPULATION_SIZE - 1);
    uniform_real_distribution<double> prob_dist(0.0, 1.0);
    vector<Individual> population(POPULATION_SIZE);
    vector<Individual> next_population;
    vector<Individual> kids_population;
    vector<Individual> combined_population;
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        population[i].values.resize(dims);
        population[i].objectives.resize(2);
        for (int j = 0; j < dims; ++j)
        {
            population[i].values[j] = init_val(rng);
        }
        // Evaluate objectives
        evaluateZDT1(population[i]);
        f_calls+=2;
        population[i].rank = 0;
        population[i].crowding_distance = 0.0;
    }
    // Calculate ranks
    paretoLayers(population);
    // Calculate crowding distances
    calculateCrowdingDistances(population);
    // --- MAIN LOOP ---
    while (f_calls < MAX_F_CALLS)
    {
        // --- 2. (crowding tournament [size=5]) SELECTION ---
        kids_population.clear();
        auto selectParent = [&]()
        {
            const Individual* best_parent = nullptr;
            for (int i = 0; i < 5; i++)
            {
                const Individual& candidate = population[pop_dist(rng)];
                if (best_parent == nullptr || candidate.rank < best_parent->rank || (candidate.rank == best_parent->rank && candidate.crowding_distance > best_parent->crowding_distance))
                    best_parent = &candidate;
            }
            return *best_parent;
        };

        // --- 3. CROSSOVER and MUTATION ---
        while (kids_population.size() < POPULATION_SIZE)
        {
            const Individual& parent1 = selectParent();
            const Individual& parent2 = selectParent();
            const Individual& parent3 = selectParent();
            Individual child;
            child.values.resize(dims);
            if (prob_dist(rng) < CROSSOVER_RATE)
            {
                //Classic intermediate real valued recombination - average from all parents
                for (int i = 0; i < dims; i++)
                {
                    child.values[i] = (parent1.values[i] + parent2.values[i] + parent3.values[i]) / 3.0;
                }
            }
            // If no crossover, just copy one parent
            else
                child.values = parent1.values;
            // Mutation - gaussian
            if (prob_dist(rng) < MUTATION_RATE)
            {
                int mutate_index = rng() % dims;
                normal_distribution<double> mutation_dist(0.0, 0.1);
                child.values[mutate_index] += mutation_dist(rng);
                if (child.values[mutate_index] < DOMAIN_MIN)
                    child.values[mutate_index] = DOMAIN_MIN;
                if (child.values[mutate_index] > DOMAIN_MAX)
                    child.values[mutate_index] = DOMAIN_MAX;
            }
            kids_population.push_back(child);
        }
        // Evaluate kids
        for (auto& kid : kids_population)
        {
            kid.objectives.resize(2);
            evaluateZDT1(kid);
            f_calls+=2;
            kid.rank = 0;
            kid.crowding_distance = 0.0;
        }
        // --- 4. FORM NEW POPULATION ---
        combined_population.clear();
        combined_population.reserve(population.size() + kids_population.size());
        combined_population.insert(combined_population.end(), population.begin(), population.end());
        combined_population.insert(combined_population.end(), kids_population.begin(), kids_population.end());
        paretoLayers(combined_population);
        calculateCrowdingDistances(combined_population);
        next_population.clear();
        map<int, vector<Individual>> fronts;
        for (const auto& ind : combined_population)
        {
            fronts[ind.rank].push_back(ind);
        }
        for (auto& [rank, front] : fronts)
        {
            // Check if we can add the whole front
            if (next_population.size() + front.size() <= POPULATION_SIZE)
                next_population.insert(next_population.end(), front.begin(), front.end());
            // if not sort by crowding distance and add the best individuals
            else
            {
                sort(front.begin(), front.end(),
                    [](const Individual& a, const Individual& b) {
                        return a.crowding_distance > b.crowding_distance;
                    });

                size_t remaining = POPULATION_SIZE - next_population.size();
                next_population.insert(next_population.end(),
                                    front.begin(),
                                    front.begin() + remaining);
                break;
            }
        }
        population = std::move(next_population);
    }
    vector<Individual> pareto_front;
    for (const auto& ind : population)
    {
        if (ind.rank == 1) {
            pareto_front.push_back(ind);
        }
    }
    return pareto_front;
}