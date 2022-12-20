#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "NumCalc.hpp"

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

std::vector<double> mulMat(std::vector<std::vector<double>> &A,
                           std::vector<double> &B) {
    int n = A.size();
    std::vector<double> ans(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            ans[i] += A[i][j] * B[j];
        }
    }
    return ans;
}

/*long long int fact(int n) {
    long long int ans = 1;
    for (int i = 2; i < n; i++) {
        ans *= i;
    }
    return ans;
}*/

double power(double b, unsigned long long e) {
    double v = 1.0;
    while(e != 0) {
        if((e & 1) != 0) {
            v *= b;
        }
        b *= b;
        e >>= 1;
    }
    return v;
}

std::vector <double> solveSys(std::vector<std::vector<double>> &sys,
                              std::vector<double> &ans) {
    inversed(sys);
    return mulMat(sys, ans);
}

std::vector <double> NeoCo(std::vector<double> points) {
    int n = points.size();

    std::vector<std::vector <double> > sys(n, std::vector <double>(n, 0));
    std::vector<double> ans(n, 0);

    long long int f = 1;
    for (int i = 0; i < n; i++) {
        if (i != 0) {
            f *= i;
        }
        for (int j = 0; j < n; j++) {
            sys[i][j] = power(points[j], i) / f;
        }
    }
    ans[1] = 1;

    return solveSys(sys, ans);
}

double to_float(std::string s) {
    int l_float = 0;
    bool flag = false, neg = false;
    double ans = 0;
    for (char i : s) {
        if (i == '-') {
            if (neg) {
                return 0;
            }
            neg = true;
            continue;
        }
        if (i == '.' || i == ',') {
            if (flag) {
                return 0;
            }
            flag = true;
        } else {
            if (i - '0' >= 0 && i - '0' <= 9) {
                ans *= 10;
                ans += i - '0';
                if (flag) {
                    l_float++;
                }
            } else {
                return 0;
            }
        }
    }
    for (int i = 0; i < l_float; i++) {
        ans /= 10;
    }
    if (neg) {
        ans *= -1;
    }
    return ans;
}

std::vector<double> strToDouble(std::string s, int n) {
    std::vector<double> ans(n, 0);
    std::vector<std::string> substring(n);
    int curr = 0;
    for (int i = 0; i < s.size(); i++) {
        if (s[i] != ' ') {
            substring[curr] += s[i];
        } else {
            curr++;
        }
    }
    for (int i = 0; i < n; i++) {
        ans[i] = to_float(substring[i]);
    }
    return ans;
}

void production() {
    std::ofstream fout("rez.txt");
    std::ifstream fin("data.txt");
    int n, N;
    std::string line;
    getline(fin, line);
    std::vector<double> data = strToDouble(line, 2);
    n = (int)data[1];
    N = (int)data[0];
    std::cout << n << ' ' << N << '\n';
    std::vector <double> points(n);
    std::vector <double> ans(n);
    for (int i = 0; i < N; i++) {
        getline(fin, line);
        points = strToDouble(line, n);
        ans = NeoCo(points);
        printVec(ans);
        for (double j : ans) {
            fout << j << ' ';
        }
        fout << '\n';
    }
    fout.close();
    fin.close();
}

/*int main() {
    single();
}*/
