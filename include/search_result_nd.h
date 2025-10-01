#pragma once
#include <Eigen/Dense>
#include <cstdlib>
#include <iostream>
#include "common.h"

enum search_method_type_nd {
    ND_BISECT,
    ND_GOLDEN_RATIO,
    ND_FIBONACCHI,
    PER_COORD_DESCEND,
    GRADIENT_DESCEND,
    ND_NONE
};

const auto search_method_string_nd = {"Bisection", "Golden ratio", "Fibonacchi", "Per coordinate descend", "Gradient Descend", "None"};

struct search_result_nd {
    search_method_type_nd type;
    Eigen::VectorXd result;
    double accuracy;
    uint64_t iterations;
    uint64_t function_probes;

    search_result_nd() : type(search_method_type_nd::ND_NONE), result(), accuracy(0.0), iterations(0), function_probes(0) {}

    search_result_nd(search_method_type_nd type, Eigen::VectorXd result, double accuracy, uint64_t iterations, uint64_t function_probes) {
        this->type = type;
        this->result = result;
        this->accuracy = accuracy;
        this->iterations = iterations;
        this->function_probes = function_probes;
    }
};

std::ostream& operator<<(std::ostream& stream, const search_result_nd& statistic);