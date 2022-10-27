#include <algorithm>
#include <exception>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "cubic_spline.hpp"

static void create_sys(std::vector<std::vector<double>>& A,
                       const std::vector<double>& a3,
                       const std::vector<double>& b3,
                       const std::vector<double>& c3);

void CubicSpline::threediag_coefs(const std::vector<double>& x,
                                  const std::vector<double>& f,
                                  std::vector<double>& a3,
                                  std::vector<double>& c3) const {
    c3[0] = (x[2] - x[1]) / ( x[2] - x[1] + x[1] - x[0] );

    for (int i = 1; i < x.size() - 1; ++i) {

        a3[i] = (x[i] - x[i-1]) / ( (x[i] - x[i-1]) + (x[i + 1] - x[i]) );
        c3[i] = (x[i+1] - x[i]) / ( x[i] - x[i-1] + x[i+1] - x[i] );
    }
}

void CubicSpline::cumlulateC(std::vector <double> &a3, std:: vector <double> &b3,
                             std::vector <double> &c3, std::vector <double> &d3) {

    const size_t M = a3.size()-1;    //M+1
    this->c[0] = 0;
    this->c[M + 1] = 0;
    std::vector <double> p(M+1, 0), q(M+1, 0);
    p[1] = -c3[0]/b3[0];
    q[1] = d3[0]/b3[0];
    for (int i = 1; i < M; i++) {
        p[i+1] = -c3[i]/(a3[i]*p[i] + b3[i]);
        q[i+1] = (d3[i] - a3[i]*q[i])/(a3[i]*p[i] + b3[i]);
    }
    c[M] = (d3[M] - a3[M]*q[M])/(p[M]*a3[M] + b3[M]);
    for (int i = M-1; i >= 1; i--) {
        this->c[i] = this->c[i+1]*p[i+1] + q[i+1];
    }
}


void CubicSpline::abd_coefs_calculation(const std::vector<double>& x,
                                        const std::vector<double>& f) {

    int N = x.size() - 1;

    for (int i = 1; i < N + 1; ++i)
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
            d3[i] = 6 * divided_difference3(x[i], x[i+1], x[i+2],
                                            f[i], f[i+1], f[i+2]);
    }
    d3[N-1] = d3[N-2];


    threediag_coefs(x, f, a3, c3);

    // Print diag matrix
    std::cout << "Coefs of diag matrix:" << std::endl;
    std::cout << "a3 : ";
    for (int i = 1; i < a3.size(); ++i)
        std::cout << a3[i] << " ";
    std::cout << std::endl;

    std::cout << "b3 : ";
    for (int i = 0; i < b3.size(); ++i)
        std::cout << b3[i] << " ";
    std::cout << std::endl;

    std::cout << "c3 : ";
    for (int i = 0; i < c3.size() - 1; ++i)
        std::cout << c3[i] << " ";
    std::cout << std::endl;
    std::cout << "###########################" << std::endl;

    std::cout << "Column:" << std::endl;
    for (int i = 0; i < d3.size(); ++i)
        std::cout << d3[i] << " ";
    std::cout << std::endl;
    std::cout << "###########################" << std::endl;

    /* cumlulateC(a3, b3, c3, d3); */

    std::vector<std::vector<double>> A(a3.size(),
                                       std::vector<double>(a3.size(), 0));
    create_sys(A, a3, b3, c3);
    std::vector<double> ans = solveSys(A, d3);

    /* for (auto i : ans) */
    /*     std::cout << i << " "; */
    /* std::cout << std::endl; */

    abd_coefs_calculation(x, f);

    this->x = std::move(x);

    // Print spline coefs
    std::cout << "Spline coefs:" << std::endl;
    for (int i = 1; i < this->a.size(); ++i)
        std::cout << this->a[i] << " ";
    std::cout << std::endl;

    for (int i = 1; i < this->a.size(); ++i)
        std::cout << this->b[i] << " ";
    std::cout << std::endl;

    for (int i = 1; i < this->a.size(); ++i)
        std::cout << this->c[i] << " ";
    std::cout << std::endl;

    for (int i = 1; i < this->a.size(); ++i)
        std::cout << this->d[i] << " ";
    std::cout << std::endl;
    std::cout << "###########################" << std::endl;
}

static void create_sys(std::vector<std::vector<double>>& A,
                       const std::vector<double>& a3,
                       const std::vector<double>& b3,
                       const std::vector<double>& c3) {

    for (auto i : A)
        i.resize(a3.size(), 0);

    for (int i = 0; i < a3.size(); ++i)
        A[i][i] = b3[i];

    for (int i = 1; i < a3.size(); ++i)
        A[i][i-1] = a3[i];

    for (int i = 0; i < a3.size() - 1; ++i)
        A[i][i+1] = c3[i];
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

std::vector <double> solveSys(std::vector <std::vector <double> > &sys,
                              std::vector <double> &ans) {

    inversed(sys);
    return mulMat(sys, ans);
}

std::vector <double> mulMat(std::vector <std::vector <double> > &A,
                            std::vector <double> &B) {
    
    int n = A.size();
    std::vector <double> ans(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            ans[i] += A[i][j] * B[j];
        }
    }
    return ans;
}

void inversed(std::vector <std::vector <double> >& A) {

    int N = A.size();
    double curr;
    double **E = new double *[N];

    for (int i = 0; i < N; i++) {
        E[i] = new double [N];
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            E[i][j] = 0.0;
            if (i == j) {
                E[i][j] = 1.0;
            }
        }
    }

    for (int k = 0; k < N; k++) {
        curr = A[k][k];
        for (int j = 0; j < N; j++) {
            A[k][j] /= curr;
            E[k][j] /= curr;
        }
        for (int i = k + 1; i < N; i++) {
            curr = A[i][k];
            for (int j = 0; j < N; j++) {
                A[i][j] -= A[k][j] * curr;
                E[i][j] -= E[k][j] * curr;
            }
        }
    }

    for (int k = N - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            curr = A[i][k];
            for (int j = 0; j < N; j++) {
                A[i][j] -= A[k][j] * curr;
                E[i][j] -= E[k][j] * curr;
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = E[i][j];
        }
    }
    for (int i = 0; i < N; i++) {
        delete [] E[i];
    }
    delete [] E;
}
