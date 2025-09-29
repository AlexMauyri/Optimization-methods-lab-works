#pragma once
#include "common.h"

enum search_method_type {
    BISECT,
    GOLDEN_RATIO,
    FIBONACCHI
};

const auto search_method_string = {"Bisection", "Golden ratio", "Fibonacchi"};

struct search_result {
    search_method_type type;
    F64 result;
    F64 accuracy;
    UI64 iterations;
    UI64 function_probes;

    search_result(search_method_type type, F64 result, F64 accuracy, UI64 iterations, UI64 function_probes) {
        this->type = type;
        this->result = result;
        this->accuracy = accuracy;
        this->iterations = iterations;
        this->function_probes = function_probes;
    }
};

std::ostream& operator<<(std::ostream& stream, const search_result& statistic);