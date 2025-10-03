#include <functional>
#include <cmath>
#include "one_dim.h"
#include "common.h"
#include "search_result.h"

search_result* bisect(const std::function<double(double)> function, double left, double right, const double eps, const uint64_t max_iterations) {
    #ifdef __DEBUG__
        std::cout << "Called one dimensional bisect method with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << "; max_iterations = " << max_iterations << '\n';
    #endif

    search_result* statistic = new search_result();
    statistic->type = search_method_type::BISECT;

    while (statistic->iterations != max_iterations && (statistic->accuracy = std::abs(right - left)) >= 2 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations + 1 << ": left = " << left
             << "; right = " << right << "; accuracy = " << statistic->accuracy << '\n';
        #endif

        statistic->result = (left + right) / 2;
        
        double x_l = statistic->result - eps / 10;
        double x_r = statistic->result + eps / 10;

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

search_result* golden_ratio(const std::function<double(double)> function, double left, double right, const double eps, const uint64_t max_iterations) {
    #ifdef __DEBUG__
        std::cout << "Called one dimensional golden ratio method with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << "; max_iterations = " << max_iterations << '\n';
    #endif

    search_result* statistic = new search_result();
    statistic->type = search_method_type::GOLDEN_RATIO;

    double x_r = left + PSI * (right - left);
    double x_l = right - PSI * (right - left);

    double y_l = function(x_l);
    double y_r = function(x_r);

    while (statistic->iterations != max_iterations && (statistic->accuracy = std::abs(right - left)) >= 2 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations + 1 << ": left = " << left
             << "; right = " << right << "; accuracy = " << statistic->accuracy << '\n';
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

search_result* fibonacchi(const std::function<double(double)> function, double left, double right, const double eps) {
    #ifdef __DEBUG__
        std::cout << "Called one dimensional fibonacchi method with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << '\n';
    #endif

    search_result* statistic = new search_result();
    statistic->type = search_method_type::FIBONACCHI;

    uint64_t fib_temp = 0, fib_1 = 1, fib_2 = 1;
    double threshold = (right - left) / eps;

    while (fib_2 < threshold) {
        fib_temp = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_temp;
        ++statistic->iterations;
    }

    double x_r = left + fib_1 * (right - left) / fib_2;
    double x_l = left + (fib_2 - fib_1) * (right - left) / fib_2;

    double y_r = function(x_r);
    double y_l = function(x_l);

    for (uint64_t iterations = statistic->iterations; iterations > 0; --iterations) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations - iterations + 1 << ": left = " << left
             << "; right = " << right << "; accuracy = " << std::abs(right - left) / 2 << '\n';
        #endif

        fib_temp = fib_2 - fib_1;
        fib_2 = fib_1;
        fib_1 = fib_temp;

        if (y_l > y_r) {
            left = x_l;
            x_l = x_r;
            y_l = y_r;
            x_r = left + fib_1 * (right - left) / fib_2;
            y_r = function(x_r);
        } else {
            right = x_r;
            x_r = x_l;
            y_r = y_l;
            x_l = left + (fib_2 - fib_1) * (right - left) / fib_2;
            y_l = function(x_l);
        }
    }

    statistic->result = (left + right) / 2;
    statistic->accuracy = std::abs(right - left) / 2;
    statistic->function_probes = statistic->iterations + 2;

    return statistic;
}