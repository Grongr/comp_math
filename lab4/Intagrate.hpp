#ifndef INTAGRATE_INCLUDED
#define INTAGRATE_INCLUDED

#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

template <typename T>
void printVec(std::vector <T> arr, std::ostream &out = std::cout, std::string sep = " ") {
    for (int i = 0; i < arr.size(); i++) {
        out << arr[i] << sep;
    }
    out << '\n';
}

double movePoint(double a0, double b0, double a1, double b1, double x);

double integral(double a, double b);

double integrateOneSeg(double a, double b, const std::function<double(double)> & func);

double integrate(double a, double b, unsigned n, const std::function<double(double)> & func);

#endif // INTAGRATE_INCLUDED