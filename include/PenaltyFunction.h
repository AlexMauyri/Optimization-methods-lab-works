#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>

#include "search_result_nd.h"

using function_nd = std::function<double(const Eigen::VectorXd&)>;

class PenaltyFunction {
protected:
    function_nd target_function;
    std::vector<function_nd> inequalities;

    virtual double target_function_with_penalty(const Eigen::VectorXd& start) const = 0;
public:
    PenaltyFunction(function_nd);

    size_t amount_of_inequalities() const;

    search_result_nd compute_minimum(const Eigen::VectorXd&) const;

    void add_inequality(function_nd);
    void delete_inequality(int);
};