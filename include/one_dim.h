#pragma once
#include "common.h"
#include "search_result.h"

search_result bisect(const std::function<double(double)> function, double left, double right, const double eps=ACCURACY, const uint64_t max_iterations=ITERS_MAX);
search_result golden_ratio(const std::function<double(double)> function, double left, double right, const double eps=ACCURACY, const uint64_t max_iterations=ITERS_MAX);
search_result fibonacchi(const std::function<double(double)> function, double left, double right, const double eps=ACCURACY);