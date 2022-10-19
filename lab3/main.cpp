#include "cubic_spline.hpp"
#include <iostream>
#include <vector>
#include <cmath>

int main() {

    std::vector<double> x = {0, 1, 2, 3};
    std::vector<double> f(x.size());

    for (int i = 0; i < f.size(); ++i)
        f[i] = std::cos(x[i]);

    CubicSpline cs(x, f);

    const double pi = std::atan(1)*4;
    double result = cs.interpolate(pi / 2);

    std::cout << result << std::endl;

    return 0;
}
