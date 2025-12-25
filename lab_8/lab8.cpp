#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;

struct Point {
    vector<double> values;
    bool dominated = false;
};

bool dominates(const Point& a, const Point& b) {
    bool strictlyBetter = false;
    for (size_t i = 0; i < a.values.size(); i++) {
        if (a.values[i] < b.values[i])
            return false;
        if (a.values[i] > b.values[i])
            strictlyBetter = true;
    }
    return strictlyBetter;
}

vector<Point> naivePareto(vector<Point> points) {
    for (auto& p : points) {
        p.dominated = false;
        for (auto& q : points) {
            if (&p != &q && dominates(q, p)) {
                p.dominated = true;
                break;
            }
        }
    }

    vector<Point> result;
    for (auto& p : points)
        if (!p.dominated)
            result.push_back(p);

    return result;
}

vector<Point> kungPareto(vector<Point>& points, int dim) {
    if (points.size() <= 1)
        return points;

    sort(points.begin(), points.end(),
         [dim](const Point& a, const Point& b) {
             return a.values[dim] > b.values[dim];
         });

    size_t mid = points.size() / 2;
    vector<Point> left(points.begin(), points.begin() + mid);
    vector<Point> right(points.begin() + mid, points.end());

    left = kungPareto(left, dim);
    right = kungPareto(right, dim);

    vector<Point> result = left;

    for (auto& r : right) {
        bool dom = false;
        for (auto& l : left) {
            if (dominates(l, r)) {
                dom = true;
                break;
            }
        }
        if (!dom)
            result.push_back(r);
    }

    return result;
}

vector<Point> generateRandomPoints(int count, int dim) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(0.0, 100.0);

    vector<Point> points(count);
    for (int i = 0; i < count; i++) {
        points[i].values.resize(dim);
        for (int j = 0; j < dim; j++)
            points[i].values[j] = dist(gen);
    }
    return points;
}

void saveCSV(const string& filename,
             const vector<Point>& all,
             const vector<Point>& pareto) {
    ofstream file(filename);
    int dim = all[0].values.size();

    for (int i = 0; i < dim; i++)
        file << "K" << i + 1 << ",";
    file << "Pareto\n";
    for (auto& p : all) {
        for (double v : p.values) file << v << ",";
        bool isPareto = false;
        for (auto& q : pareto)
            if (p.values == q.values) isPareto = true;
        file << (isPareto ? 1 : 0) << "\n";
    }
    file.close();
}

vector<vector<Point>> paretoLayers(vector<Point> points) {
    vector<vector<Point>> layers;

    while (!points.empty()) {
        auto front = naivePareto(points);
        layers.push_back(front);

        vector<Point> remaining;
        for (auto& p : points) {
            bool inFront = false;
            for (auto& f : front)
                if (p.values == f.values)
                    inFront = true;
            if (!inFront)
                remaining.push_back(p);
        }
        points = remaining;
    }
    return layers;
}

vector<Point> loadFromFile(const string& filename) {
    ifstream file(filename);
    vector<Point> points;
    double a, b;

    while (file >> a >> b) {
        Point p;
        p.values = {a, b};
        points.push_back(p);
    }
    return points;
}

int main() {

    auto points2D = generateRandomPoints(100, 2);

    auto start = chrono::high_resolution_clock::now();
    auto paretoNaive2D = naivePareto(points2D);
    auto end = chrono::high_resolution_clock::now();
    cout << "2D Naive: "
         << chrono::duration<double, milli>(end - start).count()
         << " ms\n";

    auto tmp2D = points2D;
    start = chrono::high_resolution_clock::now();
    auto paretoKung2D = kungPareto(tmp2D, 0);
    end = chrono::high_resolution_clock::now();
    cout << "2D Kung:  "
         << chrono::duration<double, milli>(end - start).count()
         << " ms\n";

    saveCSV("zad1_2D_naive.csv", points2D, paretoNaive2D);
    saveCSV("zad1_2D_kung.csv", points2D, paretoKung2D);

    auto points5D = generateRandomPoints(1000, 5);

    start = chrono::high_resolution_clock::now();
    auto paretoNaive5D = naivePareto(points5D);
    end = chrono::high_resolution_clock::now();
    cout << "5D Naive: "
         << chrono::duration<double, milli>(end - start).count()
         << " ms\n";

    auto tmp5D = points5D;
    start = chrono::high_resolution_clock::now();
    auto paretoKung5D = kungPareto(tmp5D, 0);
    end = chrono::high_resolution_clock::now();
    cout << "5D Kung:  "
         << chrono::duration<double, milli>(end - start).count()
         << " ms\n";

    cout << "2D: Naive = " << paretoNaive2D.size()
         << ", Kung = " << paretoKung2D.size() << endl;
    
    cout << "5D: Naive = " << paretoNaive5D.size()
         << ", Kung = " << paretoKung5D.size() << endl;


    saveCSV("zad1_5D_naive.csv", points5D, paretoNaive5D);
    saveCSV("zad1_5D_kung.csv", points5D, paretoKung5D);

    // zad 2
    auto filePoints = loadFromFile("MO-D3R.txt");
    auto fronts = paretoLayers(filePoints);

    ofstream file("zad2_fronty.csv");
    file << "K1,K2,Front\n";

    for (size_t i = 0; i < fronts.size(); i++) {
        for (auto& p : fronts[i]) {
            file << p.values[0] << ","
                << p.values[1] << ","
                << (i + 1) << "\n";
        }
    }

    file.close();

    return 0;
}
