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

enum class AggregatingFunction {SUM, MAX, MIN};
enum class InequalityValueFunction {INVERSE, LOG_NATURAL};

class InnerPenaltyFunction final : public PenaltyFunction {
private:
    static constexpr double EPSILON { 1e-10 };

    AggregatingFunction     agg_func;
    InequalityValueFunction ineq_value_func;

    double computeTargetFunctionWithPenalty(const Eigen::VectorXd& start) const override { return target_function(start) + innerPenalty(start); }

    double innerPenalty(const Eigen::VectorXd& start) const {
        if (inequalities.empty()) {
            return 0.0;
        }

        const std::vector<double> inequality_values = transformInequalities(start);

        switch (this->agg_func) {
            case AggregatingFunction::SUM: return std::accumulate(inequality_values.begin(), inequality_values.end(), 0.0);
            case AggregatingFunction::MAX: return *std::max_element(inequality_values.begin(), inequality_values.end());
            case AggregatingFunction::MIN: return *std::min_element(inequality_values.begin(), inequality_values.end());
        }
    }

    const std::vector<double> transformInequalities(const Eigen::VectorXd& start) const {
        std::vector<double> inequality_values;
        inequality_values.reserve(inequalities.size());
        
        for (const auto& inequality : inequalities) {
            inequality_values.push_back(inequalityValueTransform(inequality(start)));
        }

        return inequality_values;
    }

    double inequalityValueTransform(double value) const noexcept {
        switch (this->ineq_value_func) {
            case InequalityValueFunction::INVERSE: 
                if (std::abs(value) < EPSILON) return 1.0 / (EPSILON + value * value);
                else return 1.0 / (value * value);
            case InequalityValueFunction::LOG_NATURAL:
                if (std::abs(value) < EPSILON) return -std::log(EPSILON + value * value);
                else return -std::log(value * value);
        }
    }
public:
    explicit InnerPenaltyFunction(
        const function_nd& target_function,
        AggregatingFunction agg_func = AggregatingFunction::SUM,
        InequalityValueFunction ineq_value_func = InequalityValueFunction::INVERSE
    ) noexcept
        : PenaltyFunction(target_function)
        , agg_func(agg_func)
        , ineq_value_func(ineq_value_func) {}

    InnerPenaltyFunction(const InnerPenaltyFunction& other)
        : PenaltyFunction(other)
        , agg_func(other.getAggregatingFunction())
        , ineq_value_func(other.getInequalityValueFunction()) {}

    InnerPenaltyFunction(InnerPenaltyFunction&& other) noexcept
        : PenaltyFunction(std::move(other))
        , agg_func(std::move(other.agg_func))
        , ineq_value_func(std::move(other.ineq_value_func)) {}

    ~InnerPenaltyFunction() = default;

    InnerPenaltyFunction& operator=(const InnerPenaltyFunction& other) {
        if (this != &other) {
            this->agg_func = other.getAggregatingFunction();
            this->ineq_value_func = other.getInequalityValueFunction();
            PenaltyFunction::operator=(other);
        }
        return *this;
    }

    InnerPenaltyFunction& operator=(InnerPenaltyFunction&& other) noexcept {
        if (this != &other) {
            this->agg_func = std::move(other.agg_func);
            this->ineq_value_func = std::move(other.ineq_value_func);
            PenaltyFunction::operator=(std::move(other));
        }
        return *this;
    }

    inline AggregatingFunction getAggregatingFunction() const noexcept { return agg_func; }
    inline InequalityValueFunction getInequalityValueFunction() const noexcept { return ineq_value_func; }

    inline void setAggregatingFunction(AggregatingFunction agg_func) noexcept {this->agg_func = agg_func;}
    inline void setInequalityValueFunction(InequalityValueFunction ineq_value_func) noexcept {this->ineq_value_func = ineq_value_func;}
};