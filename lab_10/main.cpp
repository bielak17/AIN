//compile: g++ -fdiagnostics-color=always -g main.cpp NSGA.cpp -o main.exe
#include <iostream>
#include <ctime>
#include <fstream>

#include "NSGA.h"

using namespace std;

void save_CSV(const vector<Individual>& pareto_front, const string& filename)
{
    ofstream file(filename);
    file << "f1;f2\n";
    for (const auto& ind : pareto_front)
    {
        file << ind.objectives[0] << ";" << ind.objectives[1] << "\n";
    }
    file.close();
}

int main()
{
    runAllZDTs();
    return 0;
}