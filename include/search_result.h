#pragma once
#include <cstdlib>
#include <iostream>
#include "common.h"

enum search_method_type {
    BISECT,
    GOLDEN_RATIO,
    FIBONACCHI,
    NONE
};

const auto search_method_string = {"Bisection", "Golden ratio", "Fibonacchi", "None"};

struct search_result {
    search_method_type type;
    double result;
    double accuracy;
    uint64_t iterations;
    uint64_t function_probes;

    search_result() : type(search_method_type::NONE), result(0.0), accuracy(0.0), iterations(0), function_probes(0) {}

    search_result(search_method_type type, double result, double accuracy, uint64_t iterations, uint64_t function_probes) {
        this->type = type;
        this->result = result;
        this->accuracy = accuracy;
        this->iterations = iterations;
        this->function_probes = function_probes;
    }
};

std::ostream& operator<<(std::ostream& stream, const search_result& statistic);