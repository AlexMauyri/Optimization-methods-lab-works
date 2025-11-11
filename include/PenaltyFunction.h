#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>

#include "search_result_nd.h"
#include "common.h"

class PenaltyFunction {
private:
    const function_nd penalty_func = [this](const Eigen::VectorXd& x) {
        return this->target_function_with_penalty(x);
    };
protected:
    function_nd target_function;
    std::vector<function_nd> inequalities;

    virtual double target_function_with_penalty(const Eigen::VectorXd& start) const = 0;
public:
    PenaltyFunction(function_nd target_function) : target_function(std::move(target_function)) {}

    inline size_t amount_of_inequalities() const noexcept { return inequalities.size(); }

    inline search_result_nd compute_minimum(const Eigen::VectorXd& start, 
                                            const double eps=N_DIM_ACCURACY, 
                                            const uint64_t max_iterations=N_DIM_ITERS_MAX) const {
        return newton(penalty_func, start, eps, max_iterations);
    }

    inline void add_inequality(function_nd inequality) { inequalities.push_back(inequality); }

    inline void delete_inequality(size_t index) {
        if (index >= inequalities.size()) {
            throw std::out_of_range("Given index is out of range.");
        }

        inequalities.erase(inequalities.begin() + index);
    }

    inline const std::vector<function_nd>& get_inequalities() const { return inequalities; }

    inline void clear_inequalities() noexcept { inequalities.clear(); }

    virtual ~PenaltyFunction() = default;
};