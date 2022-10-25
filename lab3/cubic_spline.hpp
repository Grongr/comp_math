#ifndef CUBIC_SPLINE_IS_INCLUDED
#define CUBIC_SPLINE_IS_INCLUDED

// Чтоб работало

#include <vector>

class CubicSpline {

public:

    /*!
     *
     * Constructor, that calculates the coefficients of Cubic Spline
     * with (x, f) points
     *
     * @param [in] <x> points $x_i$
     * @param [in] <f> function values $f_i = f(x_i)$
     *
     */
    CubicSpline(const std::vector<double>& x,
                const std::vector<double>& f);

    /*!
     *
     * Function to calculate three diagonal coefficients
     * a3 = {a_i} and c3 = {c_i}. b3 = {b_i = 2}
     *
     * {b_0 c_0                  ...        }
     * {a_1 b_1 c_1                         }
     * {    a_2 b_2 c_2                     }
     * {                 ...                }
     * {             a_(n-1) b_(n-1) c_(n-1)}
     * {                       a_n     b_n  }
     *
     * @param [in]  <x>  points $x_i$
     * @param [in]  <f>  function values $f_i = f(x_i)$
     * @param [out] <a3> three diagonal coef
     * @param [out] <c3> three diagonal coef
     */
    void threediag_coefs(const std::vector<double>& x,
                         const std::vector<double>& f,
                         std::vector<double>& a3,
                         std::vector<double>& c3) const;
    
    /*!
     *
     * Function to calculate <c> coefficients of spline from three
     * diagonal coefficients
     *
     * @param [in] <a3> three diagonal coef
     * @param [in] <b3> three diagonal coef
     * @param [in] <c3> three diagonal coef
     * @param [in] <d3> divided differences
     */
    void cumlulateC(std::vector <double> &a3, std::vector <double> &b3,
                    std::vector <double> &c3, std::vector <double> &d3);

    /*!
     *
     * Function to calculate <a>, <b> and <d> coefficients of spline
     *
     * @param [in] <x> points $x_i$
     * @param [in] <f> function values $f_i = f(x_i)$
     */
    void abd_coefs_calculation(const std::vector<double>& x,
                               const std::vector<double>& f);

    [[nodiscard]] double divided_difference2(double x1, double x2,
                                             double f1, double f2);

    [[nodiscard]] double divided_difference3(double x1, double x2, double x3,
                                             double f1, double f2, double f3);

    [[nodiscard]] double interpolate(double x) const;

private:

    //! Coefficients of cubic spline functions
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;

    std::vector<double> x;
};

#endif //CUBIC_SPLINE_IS_INCLUDED
