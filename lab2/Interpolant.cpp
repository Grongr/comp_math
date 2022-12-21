#include "Interpolant.hpp"

static double aim(double x) {
    return std::cos(x);
}

void NewtonianInterpolant::div_differences() {

    for(size_t i = 0; i < dim - 1; i++) {
        double f = (doty[i] - doty[i + 1])/(dotx[i] - dotx[i + 1]);
        diff[i] = f;
    }

    diff[dim - 1] = diff[0];

    for(size_t i = 1; i < dim; i++) {

        for(size_t j = 0; j < dim - i - 1; j++) {

            diff[j] = (diff[j] - diff[j + 1]) /
                      (dotx[j] - dotx[j + i + 1]);
        }
        diff[dim - i - 1] = diff[0];
    }

    diff[0] = doty[0];
    reverse(diff.begin()+1, diff.end());
}

double NewtonianInterpolant::interpolate(double x, int n) const {

    if (n == dim-1) return diff[n];

    return diff[n] + (x-dotx[n])*interpolate(x, n+1);
}

double getDiff(NewtonianInterpolant& INTER, double len, int n) {
    double ans = 0;
    for (int i = 0; i < n; i++) {
        double x = len/(n-1)*i;
        ans = std::max(fabs(aim(x) - INTER.interpolate(x)), ans);
    }
    return ans;
}

void research() {
    std::ofstream fout("rez2.txt");
    std::vector<int> amount = {3, 4, 5, 7, 9, 10, 12, 15, 18, 20, 23, 26, 29, 33, 37, 41, 45, 50};
    std::vector<double> length = {256, 128, 64, 32, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125};
    printVec(amount, fout);
    printVec(length, fout);
    for (int am : amount) {
        for (double len : length) {
            std::vector<double> x(am, 0);
            std::vector<double> y(am, 0);
            for (int i = 0; i < am; i++) {
                x[i] = len/(am-1)*i;
                y[i] = aim(len/(am-1)*i);
            }
            NewtonianInterpolant INTER(x, y);
            fout << getDiff(INTER, len) << ' ';

            std::cout << getDiff(INTER, len) << ' ';
        }
        fout << std::endl;
        std::cout << std::endl;
    }
    fout.close();
}

void production() {
    std::ofstream fout("rez2.txt");
    std::vector<int> amount = {3, 5, 10};
    std::vector<double> length = {8, 4, 2, 1, 0.5, 0.25, 0.125};
    printVec(amount, fout);
    printVec(length, fout);
    for (int am : amount) {
        for (double len : length) {
            std::vector <double> x(am, 0);
            std::vector <double> y(am, 0);
            for (int i = 0; i < am; i++) {
                x[i] = len/(am-1)*i;
                y[i] = aim(x[i]);
            }
            NewtonianInterpolant INTER(x, y);
            fout << getDiff(INTER, len) << ' ';

            std::cout << getDiff(INTER, len) << ' ';
        }
        fout << std::endl;
        std::cout << std::endl;
    }
    fout.close();
}

void test1() {
    std::vector<double> x = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
    std::vector <double> y = {0, 0.033, 0.067, 0.1, 0.133, 0.166, 0.199, 0.231, 0.264, 0.296, 0.327};

    NewtonianInterpolant INTER(x, y);
    INTER.printDiff();
    std::cout << INTER.interpolate(0.95);
}

void test2() {
    std::vector<double> x = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
    std::vector <double> y(x.size(), 0);
    for (size_t i = 0; i < x.size(); i++) {
        y[i] = aim(x[i]);
    }
    NewtonianInterpolant INTER(x, y);
    std::cout << INTER.interpolate(0);
}

void KR() {
    std::vector<double> x = {0.799038105676658, -0.4999999999999999, -1.799038105676658};

    std::vector<double> y = {-0.624, 0.119, 0.123};

    NewtonianInterpolant INTER(x, y);
    for (auto i : INTER.diff) {
        std::cout << i << ' ';
    }
    std::cout << std::endl;
    std::cout << "f(0) = " << INTER.interpolate(0) << std::endl;
    std::cout << "f(pi/4) = " << INTER.interpolate(M_PI/4) << std::endl;
    std::cout << "f(pi/2) = " << INTER.interpolate(M_PI/2) << std::endl;
    std::cout << "f(pi) = " << INTER.interpolate(M_PI) << std::endl;
    std::cout << "f(-1) = " << INTER.interpolate(-1) << std::endl;
    std::cout << "Diff: " << getDiff(INTER, 5-(-2)) << std::endl;
}