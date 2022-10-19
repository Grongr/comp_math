#include <algorithm>
#include <exception>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "cubic_spline.hpp"

[[nodiscard]] double CubicSpline::divided_difference2(double x1, double x2,
                                         double f1, double f2) {
    return (f2 - f1) / (x2 - x1);
}

[[nodiscard]] double CubicSpline::divided_difference3(double x1, double x2, double x3,
                                         double f1, double f2, double f3) {
    return (divided_difference2(x2, x3, f2, f3) -
            divided_difference2(x1, x2, f1, f2)) / (x3 - x1);
}

void CubicSpline::threediag_coefs(const std::vector<double>& x,
                                                 const std::vector<double>& f,
                                                 std::vector<double>& a3,
                                                 std::vector<double>& c3) const {

    for (int i = 1; i < x.size() - 1; ++i) {

        a3[i] = (x[i] - x[i-1]) / ( (x[i] - x[i-1]) + (x[i + 1] - x[i]) );
        c3[i] = (x[i+1] - x[i]) / ( x[i] - x[i-1] + x[i+1] - x[i] );
    }
}

void CubicSpline::cumlulateC(std::vector <double> &a3, std:: vector <double> &b3,
                             std::vector <double> &c3,
                             std::vector <double> &d3) {
    const size_t M = a3.size()-1;    //M+1
    this->c[0] = 0;
    this->c[M] = 0;
    std::vector <double> p(M+1, 0), q(M+1, 0);
    p[1] = -c3[0]/2;
    q[1] = d3[0]/b3[0];
    for (int i = 1; i < M; i++) {
        p[i+1] = -c3[i]/(a3[i]*p[i] + b3[i]);
        q[i+1] = (d3[i] - a3[i]*q[i])/(a3[i]*p[i] + b3[i]);
    }
    c[M] = (d3[M] - a3[M]*q[M])/(p[M]*a3[M] + b3[M]);
    for (int i = M-1; i >= 0; i--) {
        this->c[i] = this->c[i+1]*p[i+1] + q[i+1];
    }
}


void CubicSpline::abd_coefs_calculation(const std::vector<double>& x,
                                        const std::vector<double>& f) {

    int N = x.size() - 1;

    for (int i = 0; i < N + 1; ++i)
        this->a[i] = f[i];

    this->b[1] = this->c[1] * (x[1] - x[0]) / 3 +
                 divided_difference2(x[0], x[1], f[0], f[1]);

    for (int i = 2; i < N + 1; ++i) {

        this->b[i] = (this->c[i] / 3 + this->c[i-1 / 6]) * (x[i] - x[i-1]) +
                     divided_difference2(x[i-1], x[i], f[i-1], f[i]);
        
    }

    this->d[1] = c[1] / (x[1] - x[0]);

    for (int i = 2; i < N + 1; ++i) {

        this->d[i] = (this->c[i] - this->c[i-1]) / (x[i] - x[i-1]);
    }
}



CubicSpline::CubicSpline(const std::vector<double>& x,
                         const std::vector<double>& f) {

    if (x.size() < 2)
        std::runtime_error("size of x less then 2\n");

    int N = x.size() - 1;

    this->a.resize(N + 1);
    this->b.resize(N + 1);
    this->d.resize(N + 1);
    this->c.resize(N + 1);
    
    // shitcode begin
    /* std::vector<double> a3(N - 1), b3(N - 1, 2), c3(N - 1), d3(N - 1); */
    // shitcode end
    std::vector<double> a3(N), b3(N, 2), c3(N), d3(N);

    for (int i = 0; i < N - 1; ++i) {
            d3[i] = divided_difference3(x[i], x[i+1], x[i+2],
                                 f[i], f[i+1], f[i+2]);
    }

    threediag_coefs(x, f, a3, c3);

    cumlulateC(a3, b3, c3, d3);

    abd_coefs_calculation(x, f);

    this->x = std::move(x);
}

[[nodiscard]] double CubicSpline::interpolate(double x) const {

    auto up = std::upper_bound(this->x.begin(), this->x.end(), x);

    const int i = std::min((int)(up - this->x.begin()), (int)this->x.size() - 1);

    double result = 0;

    const double tmp = x - this->x[i];

    result += a[i];
    result += b[i] * tmp;
    const double tmp2 = tmp * tmp;
    result += c[i] / 2 * tmp2;
    result += d[i] / 6 * tmp * tmp2;

    return result;
}
