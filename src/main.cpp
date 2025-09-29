#include <functional>
#include <iostream>

#include "one_dim.h"
#include "common.h"
#include "search_result.h"

static F64 test_func(const F64 x)
{
	return (x - 1) * (x - 5);
}


void lab1(std::function<F64(F64)> function) {
    F64 x0 = 0;
    F64 x1 = 10;

    std::cout << *bisect(function, x0, x1) << '\n';
    std::cout << *golden_ratio(function, x0, x1) << '\n';
    std::cout << *fibonacchi(function, x0, x1) << '\n';
}

int main() {
    lab1(test_func);

    return 0;
}