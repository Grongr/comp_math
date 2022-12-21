#ifndef NONLINEAR_INCLUDED
#define NONLINEAR_INCLUDED

#include <iostream>
#include <functional>
#include <vector>
#include <optional>
#include <fstream>
#include <limits>
#include <cmath>

enum MethodsId {NewtonId, BisectionId, IterationId};

class EquationSolver {

public:

    explicit EquationSolver(std::function<double(double)> func)
    : func{std::move(func)} {}

    static std::vector<double> GaussSolver(double** a,
                                           double* y,
                                           int n) noexcept;

    /*!
     * TODO: Write optimisation to pow
     * @param x
     * @param a
     * @return
     */
    //static double binpow(double x, int a) noexcept;

    std::function<double(double)> func;
};

class NonlinearSolver : public EquationSolver {

public:

    explicit NonlinearSolver(std::function<double(double)> func)
    : EquationSolver{std::move(func)} {}

    /*!
     * Solves equation with bisection method
     *
     * @param a - left border
     * @param b - right border
     * @param n - number of iterations
     * @return found root of equation
     */
    [[nodiscard]] double bisectionMethod(double a, double b, unsigned n) const noexcept;

    /*!
     * Solves equation with simple iteration method
     * @param inital - start approximation
     * @param tau    - method's parameter
     * @param n      - number of iterations
     * @return found root of equation
     */
    [[nodiscard]] double simpleIterationMethod(double inital, double tau,
                                               unsigned n) const noexcept;

    /*!
     * Solves equation with newton method
     * @param inital - start approximation
     * @param n      - number of iterations
     * @return found root of equation
     */
    [[nodiscard]] double newtonMethod(double inital, unsigned n) const noexcept;

    [[nodiscard]] double dfdx(double h, std::vector<double> C,
                              std::vector<double>& points, double x) const noexcept;

private:

    friend class EquationSolver;
};

#endif // NONLINEAR_INCLUDED