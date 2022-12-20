#ifndef NUMCALC_INCLUDED
#define NUMCALC_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

template <typename T>
void printVec(std::vector <T> arr,
              std::ostream &out = std::cout,
              std::string sep = " ") {

    for (int i = 0; i < arr.size(); i++) {
        out << arr[i] << sep;
    }

    out << '\n';
}

/*void print2Vec(std::vector<std::vector<double> > arr) {
    for (int i = 0; i < arr.size(); i++) {
        printVec(arr[i]);
    }
}*/

void inversed(std::vector <std::vector <double> >& A);

std::vector<double> mulMat(std::vector<std::vector<double>> &A,
                           std::vector<double> &B);

double power(double b, unsigned long long e);

std::vector <double> solveSys(std::vector<std::vector<double>> &sys,
                              std::vector<double> &ans);

std::vector <double> NeoCo(std::vector<double> points);

double to_float(std::string s);

std::vector<double> strToDouble(std::string s, int n);

void production();



#endif // NUMCALC_INCLUDED