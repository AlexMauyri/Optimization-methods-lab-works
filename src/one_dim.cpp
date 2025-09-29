#include "one_dim.h"
#include "common.h"
#include "numeric_utils.h"
#include "search_result.h"

search_result* bisect(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps, const I32 max_iterations) {
    #ifdef __DEBUG__
        std::cout << "Called bisect function with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << "; max_iterations = " << max_iterations << '\n';
    #endif
    
    UI64 total_iterations = 0;

    while ((total_iterations < max_iterations) && (abs(dif_f(right, left)) > 2 * eps)) {
        F64 x_c = sum_f(left, right) / 2;
        
        F64 x_l = x_c - eps / 10;
        F64 x_r = x_c + eps / 10;

        if (function(x_l) > function(x_r)) {
            left = x_l;
        } else {
            right = x_r;
        }

        ++total_iterations;
    }

    return new search_result(BISECT, sum_f(left, right) / 2, abs(dif_f(right, left)) / 2, total_iterations, 2 * total_iterations);
}

search_result* golden_ratio(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps, const I32 max_iterations) {
    #ifdef __DEBUG__
        std::cout << "Called golden ratio function with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << "; max_iterations = " << max_iterations << '\n';
    #endif

    F64 x_r = left + PSI * dif_f(right, left);
    F64 x_l = right - PSI * dif_f(right, left);

    F64 y_l = function(x_l);
    F64 y_r = function(x_r);

    UI64 total_iterations = 0;

    while ((total_iterations < max_iterations) && (abs(dif_f(right, left)) > 2 * eps)) {

        if (y_l > y_r) {
            x_l = x_r;
            y_l = y_r;
            x_r = left + PSI * dif_f(right, left);
            y_r = function(x_r);
        } else {
            x_r = x_l;
            y_r = y_l;
            x_l = right - PSI * dif_f(right, left);
            y_l = function(x_l);
        }

        ++total_iterations;
    }

    return new search_result(GOLDEN_RATIO, sum_f(left, right) / 2, abs(dif_f(right, left)) / 2, total_iterations, total_iterations + 2);
}

search_result* fibonacchi(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps) {
    #ifdef __DEBUG__
        std::cout << "Called fibonacchi function with parameters: left = " << left << "; right = " << right 
        << "; eps =" << eps << '\n';
    #endif

    I32 fib_temp = 0, fib_1 = 1, fib_2 = 1;
    UI64 total_iterations = 0;
    F64 threshold = dif_f(right, left) / eps;

    while (fib_2 < threshold) {
        fib_temp = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_temp;
        ++total_iterations;
    }

    F64 x_r = left + fib_1 / fib_2 * dif_f(right, left);
    F64 x_l = left + dif_f(fib_2, fib_1) / fib_2 * dif_f(right, left);

    F64 y_r = function(x_r);
    F64 y_l = function(x_l);

    for (UI64 iterations = total_iterations; iterations > 0; --iterations) {
        fib_temp = fib_2 - fib_1;
        fib_2 = fib_1;
        fib_1 = fib_temp;

        if (y_l > y_r) {
            left = x_l;
            x_l = x_r;
            y_l = y_r;
            x_r = left + fib_1 / fib_2 * dif_f(right, left);
            y_r = function(x_r);
        } else {
            right = x_r;
            x_r = x_l;
            y_r = y_l;
            x_l = left + dif_f(fib_2, fib_1) / fib_2 * dif_f(right, left);
            y_l = function(x_l);
        }

        return new search_result(FIBONACCHI, sum_f(left, right) / 2, abs(dif_f(right, left)) / 2, total_iterations, total_iterations + 2);
    }
}