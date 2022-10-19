#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "cubic_spline.hpp"

using namespace std;

double aim(double x) {
    return cos(x);
}

double err(CubicSpline &spline, double len = 3, int n = 1000) {
    double ans = 0;
    for (int i = 0; i < n; i++) {
        double x = len/(n-1)*i;
        ans = max(fabs(aim(x) - spline.interpolate(x)), ans);
    }
    return ans;
}

int main() {
    int n;
    cin >> n;
    vector <double> x(n), y(n);
    for (int i = 0; i < n; i++) {
        cin >> x[i];
        y[i] = aim(x[i]);
    }

    CubicSpline spline(x, y);
    cout << err(spline);
}
