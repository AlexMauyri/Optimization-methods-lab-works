#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>

#include "PenaltyFunction.h"
#include "multi_dim.h"
#include "search_result_nd.h"
#include "common.h"

enum class constraint_type {EQUALITY, INEQUALITY};

class OuterPenaltyFunction : public PenaltyFunction {
private:
    std::vector<function_nd> equalities;

    double target_function_with_penalty(const Eigen::VectorXd& start) const override { return target_function(start) + outer_penalty(start); }

    inline double outer_penalty(const Eigen::VectorXd& start) const { return equality_penalty(start) + inequality_penalty(start); }

    inline double equality_penalty(const Eigen::VectorXd& start) const {
        return compute_penalty(start, constraint_type::EQUALITY, 
            [this](double value) {
                return equality_value_transform(value);
            }
        );
    }

    inline double inequality_penalty(const Eigen::VectorXd& start) const {
        return compute_penalty(start, constraint_type::INEQUALITY, 
            [this](double value) {
                return inequality_value_transform(value);
            }
        );
    }

    inline double compute_penalty(const Eigen::VectorXd& start, constraint_type type, std::function<double(double)> transform) const {
        const std::vector<function_nd>& constraints = (type == constraint_type::EQUALITY)? this->equalities : this->inequalities;
        
        double accum = 0.0;
        for (const auto& constraint : constraints) {
            accum += transform(constraint(start));
        }
        return accum;
    }

    inline double equality_value_transform(double value) const { return value * value; }

    inline double inequality_value_transform(double value) const { 
        double not_negative = std::max(value, 0.0);
        return not_negative * not_negative;
    }

public:
    OuterPenaltyFunction(function_nd target_function) : PenaltyFunction(std::move(target_function)) {}

    inline size_t amount_of_equalities() const { return equalities.size(); }

    inline void add_equality(function_nd equality) { equalities.push_back(equality); }

    inline void delete_equality(size_t index) {
        if (index >= equalities.size()) {
            throw std::out_of_range("Given index is out of range.");
        }

        equalities.erase(equalities.begin() + index);
    }

    inline const std::vector<function_nd>& get_equalities() const { return equalities; }

    inline void clear_equalities() noexcept { equalities.clear(); }
};