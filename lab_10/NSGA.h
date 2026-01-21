#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<numeric>
#include<random>
#include<algorithm>
#include<map>
#include <fstream>
#include <limits>

using namespace std;

using namespace std;

enum class ZDTType { ZDT1, ZDT2, ZDT3,ZDT4, ZDT6 };

const int MAX_F_CALLS = 100001;
const int POPULATION_SIZE = 100;
const double CROSSOVER_RATE = 0.93;
const double MUTATION_RATE = 0.05;
const int DIMS[3] = {10,30,50};

struct Individual {
    vector<double> values;
    vector<double> objectives;
    int rank;
    double crowding_distance;
};

void evaluateZDT1(Individual& a);
void evaluateZDT2(Individual& a);
void evaluateZDT3(Individual& a);
void evaluateZDT6(Individual& a);
void runAllZDTs();
void evaluateIndividual(Individual& ind, ZDTType type);
bool dominates(const vector<double>& a, const vector<double>& b);
void paretoLayers(vector<Individual>& population);
void calculateCrowdingDistances(vector<Individual>& population);
vector<Individual> NSGA_II(mt19937 &rng, int dims, ZDTType zdtType);