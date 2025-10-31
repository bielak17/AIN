#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;

// Function to generate uniform distribution between 0 and 10 [U(0,10)]
void generate_uniform_dist(int count) {
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> distr(0.0, 10.0);
    ofstream outfile("uniform_distribution.txt");
    if (!outfile.is_open()) 
    {
        cerr << "Error: Could not open file for writing!" << endl;
        return;
    }
    cout << "Generating " << count << " uniform random numbers..." << endl;
    for (int i = 0; i < count; ++i) 
    {
        double value = distr(generator);
        outfile << value << endl;
    }
    outfile.close();
    cout << "Numbers saved to uniform_distribution.txt" << endl;
}

// Function to generate normal distribution with mean 0 and stddev 1 [N(0,1)]
void generate_normal_dist(int count) {
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distr(0.0, 1.0);
    ofstream outfile("normal_distribution.txt");
    if (!outfile.is_open()) 
    {
        cerr << "Error: Could not open file for writing!" << endl;
        return;
    }
    cout << "Generating " << count << " normal random numbers..." << endl;
    for (int i = 0; i < count; ++i) 
    {
        double value = distr(generator);
        outfile << value << endl;
    }
    outfile.close();
    cout << "Numbers saved to normal_distribution.txt" << endl;
}

//Task 3. Monte Carlo simulation to estimate the area of the circle.
//Rectangle is 6x6 from (-3,-3) to (3,3) and circle radius is 2 from the center (0,0).
void monte_carlo_circle_area(int samples, int number) {
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> distr(-3.0, 3.0);
    int inside_circle = 0;
    stringstream filename;
    filename << "monte_carlo_results_" << number << ".txt";
    string filename_str = filename.str();
    ofstream outfile(filename_str);
    if (!outfile.is_open()) 
    {
        cerr << "Error: Could not open file " << filename_str << " for writing!" << endl;
        return;
    }
    cout << "Running Monte Carlo simulation with " << samples << " samples..." << endl;
    for (int i = 0; i < samples; ++i) {
        double x = distr(generator);
        double y = distr(generator);
        bool inside = (x * x + y * y <= 4.0);
        if (inside) 
        {
            inside_circle++;
        }
        outfile << x << ";" << y << ";" << (inside ? 1 : 0) << endl;
    }
    outfile.close();
    cout << "Results saved to " << filename_str << endl;
    double rectangle_area = 36.0;
    double circle_area_estimate = (static_cast<double>(inside_circle) / samples) * rectangle_area;
    cout << "Estimated area of the circle: " << circle_area_estimate << endl;
    cout << "Actual area of the circle: " << 3.14159 * 4 << endl;
    cout << "Error: " << fabs(circle_area_estimate - (3.14159 * 4)) << endl;
}

int main(){
    //Task 1 and task 2. Creating 2 files with 10,000 numbers of each distribution.
    generate_uniform_dist(10000);
    generate_normal_dist(10000);
    //Task 3. Monte Carlo simulation with 100 / 1000 / 10000 / 100000 samples.
    monte_carlo_circle_area(100,1);
    monte_carlo_circle_area(1000,2);
    monte_carlo_circle_area(10000,3);
    monte_carlo_circle_area(100000,4);
    return 0;
}