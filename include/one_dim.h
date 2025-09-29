#pragma once
#include "common.h"
#include "search_result.h"

search_result* bisect(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps=ACCURACY, const I32 max_iterations=ITERS_MAX);
search_result* golden_ratio(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps=ACCURACY, const I32 max_iterations=ITERS_MAX);
search_result* fibonacchi(std::function<F64(F64)> function, F64 left, F64 right, const F64 eps=ACCURACY);