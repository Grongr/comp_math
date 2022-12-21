#include <iostream>
#include "NonlinearEquations.hpp"

/*!
 *  Test function
 * @param x
 * @return
 */
double f(double x) {
    return x*x - 4;
}

/*!
 *  Kepler equation
 * @param x
 * @return f(x)
 */
double Kepler(double x) {
    return x - 0.1 * std::sin(x) - 3.14/4;
}

/*!
 * Second equation
 * @param x
 * @return f(x)
 */
double eq2(double x) {
    return std::tan(x) - 4*x/3.14;
}

/*!
 * Third equation
 * @param x
 * @return f(x)
 */
double eq3(double x) {
    return log(std::cosh(x));
}

static void solve_the_equation(std::string file_x, std::string file_y,
                               NonlinearSolver& slv, MethodsId id,
                               double accur) noexcept {

    std::ofstream filex, filey;
    filex.open(file_x);
    filey.open(file_y);

    if (id == NewtonId) {

        for (int i = 0; i < 10; ++i) {
            double ans = slv.newtonMethod(0, i);
            filex << i << std::endl;
            filey << std::abs(accur - ans) << std::endl;
        }

    } else if (id == BisectionId) {

        for (int i = 0; i < 10; ++i) {
            double ans = slv.bisectionMethod(0, 1, i + 1);
            filex << i + 1 << std::endl;
            filey << std::abs(accur - ans);
        }
    } else if (id == IterationId) {

        for (int i = 0; i < 10; ++i) {
            double ans = slv.simpleIterationMethod(1,
                          -0.5 / double(i + 1), 10);
            filex << -0.5 / double(i + 1) << std::endl;
            filey << std::abs(accur - ans) << std::endl;
        }
    }

    filex.close();
    filey.close();
}


int main() {

    NonlinearSolver slvK(Kepler);
    NonlinearSolver slv2(eq2);
    NonlinearSolver slv3(eq3);

    /* Kepler Solutions BEGIN */
    solve_the_equation("Kepler_Newton_x.dat",
                       "Kepler_Newton_y.dat",
                       slvK, NewtonId, 0.861265);

    solve_the_equation("Kepler_Bin_x.dat",
                       "Kepler_Bin_y.dat",
                       slvK, BisectionId, 0.81265);

    solve_the_equation("Kepler_Iter_x.dat",
                       "Kepler_Iter_y.dat",
                       slvK, IterationId, 0.81265);
    /* Kepler Solutions END */

    /* Func 2 Solutions BEGIN */
    solve_the_equation("Tg_Newton_x.dat",
                       "Tg_Newton_y.dat",
                       slv2, NewtonId, 0.78539816);

    solve_the_equation("Tg_Bin_x.dat",
                       "Tg_Bin_y.dat",
                       slv2, BisectionId, 0.78539816);

    solve_the_equation("Tg_Iter_x.dat",
                       "Tg_Iter_y.dat",
                       slv2, IterationId, 0.78539816);
    /* Func 2 Solutions END */

    /* Func 3 Solutions BEGIN */
    solve_the_equation("Ln_Newton_x.dat",
                       "Ln_Newton_y.dat",
                       slv3, NewtonId, 0);

    solve_the_equation("Ln_Bin_x.dat",
                       "Ln_Bin_y.dat",
                       slv3, BisectionId, 0);

    solve_the_equation("Ln_Iter_x.dat",
                       "Ln_Iter_y.dat",
                       slv3, IterationId, 0);
    /* Func 3 Solutions END */

    return 0;
}