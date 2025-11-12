#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>

#include "search_result_nd.h"
#include "common.h"

class PenaltyFunction {
private:
    const function_nd penalty_func = [this](const Eigen::VectorXd& x) {
        return this->computeTargetFunctionWithPenalty(x);
    };
protected:
    function_nd                 target_function;
    std::vector<function_nd>    inequalities;

    virtual double computeTargetFunctionWithPenalty(const Eigen::VectorXd& start) const = 0;
public:
    explicit PenaltyFunction(const function_nd& target_function) noexcept
        : target_function(target_function) {}

    PenaltyFunction(const PenaltyFunction& other) noexcept 
        : target_function(other.getTargetFunction())
        , inequalities(other.getInequalities()) {}

    PenaltyFunction(PenaltyFunction&& other) noexcept 
        : target_function(other.getTargetFunction())
        , inequalities(other.getInequalities()) {}

    ~PenaltyFunction() { inequalities.clear(); }
    
    PenaltyFunction& operator=(const PenaltyFunction& other) noexcept {
        if (this != &other) {
            this->target_function = other.getTargetFunction();
            this->inequalities.assign(other.getInequalities().begin(), other.getInequalities().end());
        }
        return *this;
    }

    PenaltyFunction& operator=(PenaltyFunction&& other) noexcept {
        if (this != &other) {
            this->target_function = other.getTargetFunction();
            this->inequalities.assign(other.getInequalities().begin(), other.getInequalities().end());
        }
        return *this;
    }

    search_result_nd computeMinimum(const Eigen::VectorXd& start, 
                                    const double eps = N_DIM_ACCURACY, 
                                    const uint64_t max_iterations = N_DIM_ITERS_MAX) const {
        return newton(penalty_func, start, eps, max_iterations);
    }

    inline size_t getAmountOfInequalities() const noexcept { return inequalities.size(); }

    inline function_nd getTargetFunction() const { return target_function; }

    inline const std::vector<function_nd>& getInequalities() const { return inequalities; }

    inline void setTargetFunction(function_nd target_function) { this->target_function = target_function; }

    inline void addInequality(function_nd inequality) { inequalities.push_back(inequality); }

    inline void deleteInequality(size_t index) {
        if (index >= inequalities.size()) {
            throw std::out_of_range("Given index is out of range.");
        }

        inequalities.erase(inequalities.begin() + index);
    }

    inline void clearInequalities() noexcept { inequalities.clear(); }
};