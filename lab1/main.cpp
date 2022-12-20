#include <iostream>
#include <vector>

#include "NumCalc.hpp"

static void single() {

    int n;
    std::cin >> n;
    std::vector<double> points(n);
    for (int i = 0; i < n; i++) {
        std::cin >> points[i];
    }
    printVec(NeoCo(points), std::cout, ", ");
}

int main() {

    single();

    return 0;
}