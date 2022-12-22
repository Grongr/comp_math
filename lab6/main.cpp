#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iomanip>

using Time = double;

#define GACHI_R 40e6

double mu = 3.9860044158e14;

struct State {

    std::vector<double> state;
    Time t;

    State(int dim) { state.resize(dim, 0); }
};

[[nodiscard]] double nvec_norm(const std::vector<double>& vec) {

    double sum = 0;

    for (size_t i = 0; i < vec.size(); ++i)
        sum += vec[i] * vec[i];

    return std::sqrt(sum);
}

[[nodiscard]] std::vector<double> rightSide (const Time time, const std::vector<double>& state) noexcept {

    std::vector<double> dstate(6, 0);
    std::vector<double> r(3);

    for (size_t i = 0; i < r.size(); ++i)
        r[i] = state[i];
    
    double norm3 = nvec_norm(r);
    norm3 *= norm3 * norm3;


    std::vector<double> tmp(r);

    for (size_t i = 0; i < tmp.size(); ++i)
        tmp[i] *= -mu / norm3;

    for (size_t i = 0; i < 3; ++i) {

        dstate[i] = state[3 + i];
        dstate[3 + i] = tmp[i];
    }

    return dstate;
}

std::vector<double> operator*(const double value, const std::vector<double>& vec) {

    std::vector<double> result(vec);

    for (size_t i = 0; i < result.size(); ++i) {

        result[i] *= value;
    }

    return result;
}

std::vector<double> operator*(const std::vector<double>& vec, const double value) {
    return value * vec;
}

std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double> vec2) {

    if (vec1.size() != vec2.size())
        throw std::runtime_error("Vectors must be the same length");

    std::vector<double> result(vec1);

    for (size_t i = 0; i < result.size(); ++i)
        result[i] += vec2[i];

    return result;
}



double implicitRk(const State& initial, Time step, unsigned iteration) {
    
    State state_k = initial;

    double r_res = 0;
    for (unsigned i = 0; i < iteration; ++i) {

        auto k1 = rightSide(state_k.t,            state_k.state);
        auto k2 = rightSide(state_k.t + step / 2, state_k.state + step/2 * k1);
        auto k3 = rightSide(state_k.t + step / 2, state_k.state + step/2 * k2);
        auto k4 = rightSide(state_k.t + step,     state_k.state + step * k3);

        state_k.t += step;
        state_k.state = state_k.state + step/6 * (k1 + 2*k2 + 2*k3 + k4);

        /* for (size_t i = 0; i < state_k.state.size(); ++i) */
        /*     std::cout << state_k.state[i] << " "; */

        /* std::cout << state_k.state[0] << "," << state_k.state[1] << std::endl; */

        double tmp = std::sqrt(state_k.state[0] * state_k.state[0] + 
                               state_k.state[1] * state_k.state[1]);

        tmp = std::abs(tmp - GACHI_R);

        r_res = tmp ? tmp > r_res : r_res;
    }

    return r_res;
}

int main() {

    State st(6);
    
    st.state[0] = 0;
    st.state[1] = 40e6;
    st.state[2] = 0;
    st.state[3] = 0;
    st.state[4] = 0;
    st.state[5] = 0;

    double norm = nvec_norm(st.state);
    st.state[3] = std::sqrt(mu / norm);

    Time step = 1;

    std::vector<double> r_res(30);
    for (step = 100; step <= 1000; step += 100) {

        unsigned iteration = 100000 / step;

        double res = implicitRk(st, step, iteration);
        r_res[(int)step/100] = res;
    }

    for (int i = 1; i < 10; ++i) {

        std::cout << std::setprecision(30);
        std::fixed(std::cout);
        std::cout << i*100 << "," << r_res[i]/GACHI_R << std::endl;
    }

    return 0;
}
