#ifndef INTERPOLANT_INCLUDED
#define INTERPOLANT_INCLUDED

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

template <typename T>
void printVec(std::vector<T> arr,
              std::ostream &out = std::cout) {

    for (size_t i = 0; i < arr.size(); i++) {
        out << arr[i] << ' ';
    }
    out << std::endl;
}

class NewtonianInterpolant {
public:
    std::vector<double> dotx, doty;
    std::vector<double> diff;
    size_t dim;

    void div_differences();

    void printDiff() const { printVec(diff); }

    double interpolate(double x, int n = 0) const;

    NewtonianInterpolant(std::vector<double> &x, std::vector<double> &y)
                        : dotx{x}, doty{y} {
        dim = x.size();
        diff.resize(dim, 0);

        div_differences();
    }
};

double getDiff(NewtonianInterpolant &INTER, double len, int n = 1000000);

void research();

void production();

void test1();

void test2();

void KR();

#endif // INTERPOLANT_INCLUDED