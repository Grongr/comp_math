#include "NonlinearEquations.hpp"

#include <cmath>

[[nodiscard]] static std::vector<double> firstDer(const std::vector<double>& points) noexcept {

    int n = points.size();
    double** A = new double* [n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        for (int j = 0; j < n; j++) {
            A[i][j] = std::pow(points[j], i);
        }
    }

    double* y = new double[n];
    for (int i = 0; i < n; i++) {
        y[i] = 0;
    };
    y[1] = 1;

    std::vector<double> x = EquationSolver::GaussSolver(A, y, n);

    delete[] y;

    for (int i = 0; i < n; i++)
        delete[] A[i];
    delete[] A;

    return x;
}

std::vector<double> EquationSolver::GaussSolver(double** a,
                                                double* y,
                                                int n) noexcept {

    double max;
    int k, index;
    const double eps = std::numeric_limits<double>::min();
    std::vector<double> x;
    k = 0;
    while (k < n)
    {
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(a[i][k]) > max)
            {
                max = abs(a[i][k]);
                index = i;
            }
        }

        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }

        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;

        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (abs(temp) < eps) continue;
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue;
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }

    for (k = n - 1; k >= 0; k--)
    {
        x.insert(x.begin(), y[k]);
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * y[k];
    }
    return x;
}

double NonlinearSolver::dfdx(double h, std::vector<double> C,
                             std::vector<double>& points, double x) const noexcept {
    double y = 0;
    for (int i = 0; i < points.size(); i++) {
        y += this->func(x + points[i] * h) * C[i] / h;
    }
    return y;
}

[[nodiscard]] double NonlinearSolver::bisectionMethod(double a, double b,
                                     unsigned n) const noexcept {
    double c = 0; // корень
    for (int i = 0; i < n; i++) {
        c = (a + b) / 2;
        if (this->func(a) * this->func(c) > 0) a = c;
        else                                   b = c;
    }
    return c;
}

[[nodiscard]] double NonlinearSolver::simpleIterationMethod(double inital, double tau,
                                                            unsigned n) const noexcept {
    double  x = inital;
    for (int i = 0; i < n; i++) {
        x = x + tau * this->func(x);
    }
    return x;
}

[[nodiscard]] double NonlinearSolver::newtonMethod(double inital, unsigned n) const noexcept {
    std::vector<double> points = {-2, -1, 0, 1, 2};
    std::vector<double> c = firstDer(points);
    double x = inital;
    for (int i = 0; i < n; i++) {
        x = x - (this->func(x) / dfdx(0.001, c, points, x));
    }
    return x;
}
