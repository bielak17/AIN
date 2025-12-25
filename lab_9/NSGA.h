#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<numeric>
#include<random>
#include<algorithm>
#include<map>

using namespace std;

const int MAX_F_CALLS = 20000;
const double DOMAIN_MIN = 0.0;
const double DOMAIN_MAX = 1.0;
const int POPULATION_SIZE = 100;
const int MAX_GENERATIONS = 100;
const double CROSSOVER_RATE = 0.9;
const double MUTATION_RATE = 0.03;

struct Individual {
    vector<double> values;
    vector<double> objectives;
    int rank;
    double crowding_distance;
};

void evaluateZDT1(Individual& a);
bool dominates(const vector<double>& a, const vector<double>& b);
void paretoLayers(vector<Individual>& population);
void calculateCrowdingDistances(vector<Individual>& population);
vector<Individual> NSGA_II(std::mt19937 &rng, int dims=2);