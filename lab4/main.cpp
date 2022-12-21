#include "Intagrate.hpp"

double func(double n) {
    return std::cos(n);
}

int main() {
    std::function<double(double)> f = func;
    std::ofstream fout("Integrate_err.txt");
    double a=0, b=10, INTEGRAL;
    int N = 10000;
    //std::cin >> a >> b;
    std::cout.precision(15);
    INTEGRAL = integral(a, b);
    for (int i = 1; i < N+1 ; i++) {
        fout << fabs(integrate(a, b, i, f) - INTEGRAL) << ' ';
        if (i % 100 == 0) {
            std::cout << i << "/" << N << '\n';
        }
    }
    fout.close();
}