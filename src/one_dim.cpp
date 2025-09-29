#include "one_dim.h"
#include "common.h"
#include "numeric_utils.h"
#include "search_result.h"

search_result* bisect(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps, const UI64 max_iterations) {
    #ifdef __DEBUG__
        std::cout << "Called bisect method with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << "; max_iterations = " << max_iterations << '\n';
    #endif

    search_result* statistic = new search_result();
    statistic->type = search_method_type::BISECT;

    while (statistic->iterations != max_iterations && (statistic->accuracy = std::abs(right - left)) >= 2 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations + 1 << ": left = " << left << "; right = " << right << '\n';
        #endif

        statistic->result = (left + right) / 2;
        
        F64 x_l = statistic->result - eps / 10;
        F64 x_r = statistic->result + eps / 10;

        if (function(x_l) > function(x_r)) {
            left = x_l;
        } else {
            right = x_r;
        }

        ++statistic->iterations;
    }

    statistic->function_probes = 2 * statistic->iterations;

    return statistic;
}

search_result* golden_ratio(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps, const UI64 max_iterations) {
    #ifdef __DEBUG__
        std::cout << "Called golden ratio method with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << "; max_iterations = " << max_iterations << '\n';
    #endif

    search_result* statistic = new search_result();
    statistic->type = search_method_type::GOLDEN_RATIO;

    F64 x_r = left + PSI * (right - left);
    F64 x_l = right - PSI * (right - left);

    F64 y_l = function(x_l);
    F64 y_r = function(x_r);

    while (statistic->iterations != max_iterations && (statistic->accuracy = std::abs(right - left)) >= 2 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations + 1 << ": left = " << left << "; right = " << right << '\n';
        #endif
        if (y_l > y_r) {
            left = x_l;
            x_l = x_r;
            y_l = y_r;
            x_r = left + PSI * (right - left);
            y_r = function(x_r);
        } else {
            right = x_r;
            x_r = x_l;
            y_r = y_l;
            x_l = right - PSI * (right - left);
            y_l = function(x_l);
        }

        ++statistic->iterations;
    }

    statistic->function_probes = statistic->iterations + 2;
    statistic->result = (left + right) / 2;

    return statistic;
}

search_result* fibonacchi(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps) {
    #ifdef __DEBUG__
        std::cout << "Called fibonacchi method with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << '\n';
    #endif

    search_result* statistic = new search_result();
    statistic->type = search_method_type::FIBONACCHI;

    I32 fib_temp = 0, fib_1 = 1, fib_2 = 1;
    F64 threshold = (right - left) / eps;

    while (fib_2 < threshold) {
        fib_temp = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_temp;
        ++statistic->iterations;
    }

    F64 x_r = left + static_cast<double>(fib_1) / fib_2 * (right - left);
    F64 x_l = left + static_cast<double>(fib_2 - fib_1) / fib_2 * (right - left);

    F64 y_r = function(x_r);
    F64 y_l = function(x_l);

    for (UI64 iterations = statistic->iterations; iterations > 0; --iterations) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations - iterations + 1 << ": left = " << left << "; right = " << right << '\n';
        #endif

        fib_temp = fib_2 - fib_1;
        fib_2 = fib_1;
        fib_1 = fib_temp;

        if (y_l > y_r) {
            left = x_l;
            x_l = x_r;
            y_l = y_r;
            x_r = left + static_cast<double>(fib_1) / fib_2 * (right - left);
            y_r = function(x_r);
        } else {
            right = x_r;
            x_r = x_l;
            y_r = y_l;
            x_l = left + static_cast<double>(fib_2 - fib_1) / fib_2 * (right - left);
            y_l = function(x_l);
        }
    }

    statistic->result = (left + right) / 2;
    statistic->accuracy = std::abs(right - left) / 2;
    statistic->function_probes = statistic->iterations + 2;

    return statistic;
}