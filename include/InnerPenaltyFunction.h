#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

#include "PenaltyFunction.h"
#include "multi_dim.h"
#include "search_result_nd.h"
#include "common.h"

enum class aggregating_function {SUM, MAX, MIN};
enum class inequality_value_function {INVERSE, LOG_NATURAL};

class InnerPenaltyFunction : public PenaltyFunction {
private:
    double target_function_with_penalty(const Eigen::VectorXd& start) const override { return target_function(start) + inner_penalty(start); }

    double inequality_value_transform(double value) const {
        constexpr double epsilon = 1e-10;

        switch (this->ineq_value_func) {
            case inequality_value_function::INVERSE: 
                if (std::abs(value) < epsilon) return std::copysign(1.0 / epsilon, value);
                else return 1.0 / value;
            case inequality_value_function::LOG_NATURAL:
                if (value < epsilon) return std::log(epsilon);
                else return std::log(value);
            default: throw std::invalid_argument("Unknown type of inequality value function.");
        }
    }

    double inner_penalty(const Eigen::VectorXd& start) const {
        if (inequalities.empty()) {
            return 0.0;
        }

        std::vector<double> inequality_values;
        inequality_values.reserve(inequalities.size());
        
        for (const auto& inequality : inequalities) {
            inequality_values.push_back(inequality_value_transform(inequality(start)));
        }

        switch (this->agg_func) {
            case aggregating_function::SUM: return std::accumulate(inequality_values.begin(), inequality_values.end(), 0.0);
            case aggregating_function::MAX: return *std::max_element(inequality_values.begin(), inequality_values.end());
            case aggregating_function::MIN: return *std::min_element(inequality_values.begin(), inequality_values.end());
            default: throw std::invalid_argument("Unknown type of aggregating function.");
        }
    }

    aggregating_function agg_func;
    inequality_value_function ineq_value_func;
public:
    InnerPenaltyFunction(function_nd target_function, 
                        aggregating_function agg_func = aggregating_function::SUM, 
                        inequality_value_function ineq_value_func = inequality_value_function::INVERSE) 
                        : PenaltyFunction(std::move(target_function))
                        , agg_func(agg_func)
                        , ineq_value_func(ineq_value_func) {}

    inline void set_aggregating_function(aggregating_function agg_func) noexcept {this->agg_func = agg_func;}
    inline void set_inequality_value_function(inequality_value_function ineq_value_func) noexcept {this->ineq_value_func = ineq_value_func;}

    inline aggregating_function get_aggregating_function() const noexcept { return agg_func; }
    inline inequality_value_function get_inequality_value_function() const noexcept { return ineq_value_func; }
};