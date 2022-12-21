#include "Intagrate.hpp"

std::vector <double> P = {-sqrt((double)3/(double)7 + (double)2/(double)7*sqrt((double)6/(double)5)),
                          -sqrt((double)3/(double)7 - (double)2/(double)7*sqrt((double)6/(double)5)),
                          sqrt((double)3/(double)7 - (double)2/(double)7*sqrt((double)6/(double)5)),
                          sqrt((double)3/(double)7 + (double)2/(double)7*sqrt((double)6/(double)5))};
std::vector <double> W = {(18-sqrt(30))/36, (18+sqrt(30))/36 , (18+sqrt(30))/36, (18-sqrt(30))/36};

double movePoint(double a0, double b0, double a1, double b1, double x) {
    return a1 + (x-a0)*(b1-a1)/(b0-a0);
}

double integral(double a, double b) {
    return sin(b) - sin(a);
}

double integrateOneSeg(double a, double b, const std::function<double(double)> & func) {
    double ans = 0;
    double bma2 = (b-a)/2, apb2 = (a+b)/2;
    for (int i = 0; i < 4; i++) {
        ans += func(bma2*P[i] + apb2)*W[i];
    }
    return ans*(b-a)/2;
}

double integrate(double a, double b, unsigned n, const std::function<double(double)> & func) {
    double ans = 0, l = (b-a)/n;
    for (unsigned i = 0; i < n; i++) {
        ans += integrateOneSeg(a+l*i, a+l*(i+1), func);
    }
    return ans;
}

