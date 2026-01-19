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

void evaluateZDT2(Individual& a) 
{ 
    double f1 = a.values[0]; 
    double g = 1+9*a.values[1];
    double h = 1-(f1/g)*(f1/g); 
    double f2 = g*h; 
    a.objectives[0] = f1; 
    a.objectives[1] = f2; 
} 
void evaluateZDT3(Individual& a) 
{ 
    double f1 = a.values[0];
    double g = 1+9*a.values[1];
    double h = 1-sqrt(f1/g)-(f1/g)*sin(10*M_PI*f1);
    double f2 = g*h;
    a.objectives[0] = f1; 
    a.objectives[1] = f2; 
} 

void evaluateZDT4(Individual& a)
{
    double f1 = a.values[0];
    double g = 1.0 + 10.0 * (a.values.size() - 1);

    for (size_t i = 1; i < a.values.size(); i++)
        g += a.values[i] * a.values[i] - 10.0 * cos(4.0 * M_PI * a.values[i]);

    double h = 1.0 - sqrt(f1 / g);
    a.objectives[0] = f1;
    a.objectives[1] = g * h;
}

void evaluateZDT6(Individual& a) 
{ 
    double f1 = 1 - exp(-4*a.values[0])*pow(sin(6*M_PI*a.values[0]),6);
    double g = 1+9*sqrt(sqrt(a.values[1])); 
    double h = 1-(f1/g)*(f1/g); 
    double f2 = g*h; 
    a.objectives[0] = f1; 
    a.objectives[1] = f2; 
}

void evaluateIndividual(Individual& ind, ZDTType type)
{
    switch (type)
    {
        case ZDTType::ZDT1: evaluateZDT1(ind); break;
        case ZDTType::ZDT2: evaluateZDT2(ind); break;
        case ZDTType::ZDT3: evaluateZDT3(ind); break;
        case ZDTType::ZDT4: evaluateZDT4(ind); break;
        case ZDTType::ZDT6: evaluateZDT6(ind); break;
    }
}

string zdtName(ZDTType type)
{
    switch (type)
    {
        case ZDTType::ZDT1: return "ZDT1";
        case ZDTType::ZDT2: return "ZDT2";
        case ZDTType::ZDT3: return "ZDT3";
        case ZDTType::ZDT4: return "ZDT4";
        case ZDTType::ZDT6: return "ZDT6";
    }
    return "";
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

        for (auto* p : front)
            p->rank = currentRank;

        vector<Individual*> newRemaining;
        for (auto* p : remaining)
        {
            if (find(front.begin(), front.end(), p) == front.end())
                newRemaining.push_back(p);
        }

        remaining = newRemaining;
        currentRank++;
    }
}

void runAllZDTs()
{
    random_device rd;
    mt19937 rng(rd());

    vector<ZDTType> zdt_functions = {
        ZDTType::ZDT1,
        ZDTType::ZDT2,
        ZDTType::ZDT3,
        ZDTType::ZDT4,
        ZDTType::ZDT6
    };

    for (ZDTType zdt : zdt_functions)
    {
        for (int dims : DIMS)
        {
            NSGA_II(rng, dims, zdt);
        }
    }
}

void calculateCrowdingDistances(vector<Individual>& population)
{
    for (auto& ind : population)
        ind.crowding_distance = 0.0;

    int maxRank = 0;
    for (const auto& ind : population)
        maxRank = max(maxRank, ind.rank);

    for (int rank = 1; rank <= maxRank; rank++)
    {
        vector<Individual*> front;
        for (auto& ind : population)
            if (ind.rank == rank)
                front.push_back(&ind);

        for (int m = 0; m < 2; m++)
        {
            sort(front.begin(), front.end(),
                [m](const Individual* a, const Individual* b)
                {
                    return a->objectives[m] < b->objectives[m];
                });

            front.front()->crowding_distance =
            front.back()->crowding_distance =
                numeric_limits<double>::infinity();

            double fmin = front.front()->objectives[m];
            double fmax = front.back()->objectives[m];

            if (fmax - fmin == 0.0)
                continue;

            for (size_t i = 1; i < front.size() - 1; i++)
            {
                front[i]->crowding_distance +=
                    (front[i + 1]->objectives[m] -
                     front[i - 1]->objectives[m]) / (fmax - fmin);
            }
        }
    }
}



void saveParetoFront(const vector<Individual>& population,
                     int iteration,
                     int dims,
                     ZDTType type)
{
    string filename =
        "pareto_iter_" + to_string(iteration) + "_" +
        zdtName(type) + "_dims_" + to_string(dims) + ".txt";

    ofstream file(filename);
    for (const auto& ind : population)
        if (ind.rank == 1)
            file << ind.objectives[0] << " "
                 << ind.objectives[1] << "\n";
}


vector<Individual> NSGA_II(mt19937& rng, int dims, ZDTType zdtType)
{
    int f_calls = 0;
    int iteration = 0;

    uniform_real_distribution<double> prob_dist(0.0, 1.0);
    uniform_int_distribution<int> pop_dist(0, POPULATION_SIZE - 1);

    vector<Individual> population(POPULATION_SIZE);
    vector<Individual> kids_population, combined_population, next_population;

    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        population[i].values.resize(dims);
        population[i].objectives.resize(2);

        for (int j = 0; j < dims; j++)
        {
            if (zdtType == ZDTType::ZDT4 && j > 0)
                population[i].values[j] =
                    uniform_real_distribution<double>(-5.0, 5.0)(rng);
            else
                population[i].values[j] =
                    uniform_real_distribution<double>(0.0, 1.0)(rng);
        }

        evaluateIndividual(population[i], zdtType);
        f_calls += 2;
    }

    paretoLayers(population);
    calculateCrowdingDistances(population);


    while (f_calls < MAX_F_CALLS)
    {
        iteration++;
        kids_population.clear();

        auto selectParent = [&]()
        {
            const Individual* best = nullptr;
            for (int i = 0; i < 5; i++)
            {
                const Individual& cand = population[pop_dist(rng)];
                if (!best || cand.rank < best->rank ||
                    (cand.rank == best->rank &&
                     cand.crowding_distance > best->crowding_distance))
                    best = &cand;
            }
            return *best;
        };

        while (kids_population.size() < POPULATION_SIZE)
        {
            const Individual& p1 = selectParent();
            const Individual& p2 = selectParent();
            const Individual& p3 = selectParent();

            Individual child;
            child.values.resize(dims);

            if (prob_dist(rng) < CROSSOVER_RATE)
            {
                for (int i = 0; i < dims; i++)
                    child.values[i] =
                        (p1.values[i] + p2.values[i] + p3.values[i]) / 3.0;
            }
            else
                child.values = p1.values;

            if (prob_dist(rng) < MUTATION_RATE)
            {
                int idx = rng() % dims;
                normal_distribution<double> mut(0.0, 0.1);
                child.values[idx] += mut(rng);

                if (zdtType == ZDTType::ZDT4 && idx > 0)
                {
                    child.values[idx] =
                        max(-5.0, min(5.0, child.values[idx]));
                }
                else
                {
                    child.values[idx] =
                        max(0.0, min(1.0, child.values[idx]));
                }
            }

            child.objectives.resize(2);
            evaluateIndividual(child, zdtType);
            f_calls += 2;
            kids_population.push_back(child);
        }

        combined_population = population;
        combined_population.insert(combined_population.end(),
                                   kids_population.begin(),
                                   kids_population.end());

        paretoLayers(combined_population);
        calculateCrowdingDistances(combined_population);

        next_population.clear();
        map<int, vector<Individual>> fronts;
        for (const auto& ind : combined_population)
            fronts[ind.rank].push_back(ind);

        for (auto& [rank, front] : fronts)
        {
            if (next_population.size() + front.size() <= POPULATION_SIZE)
                next_population.insert(next_population.end(),
                                       front.begin(), front.end());
            else
            {
                sort(front.begin(), front.end(),
                     [](const Individual& a, const Individual& b)
                     {
                         return a.crowding_distance >
                                b.crowding_distance;
                     });

                size_t remain =
                    POPULATION_SIZE - next_population.size();
                next_population.insert(next_population.end(),
                                       front.begin(),
                                       front.begin() + remain);
                break;
            }
        }

        population = next_population;

        if (iteration == 20 || iteration == 50 ||
            iteration == 100 || iteration == 500)
        {
            saveParetoFront(population, iteration, dims, zdtType);
        }
    }

    vector<Individual> pareto_front;
    for (const auto& ind : population)
        if (ind.rank == 1)
            pareto_front.push_back(ind);

    return pareto_front;
}


