#pragma once

#include <Eigen/Dense>
#include "common.h"

class SimplexResult final {
private:
    const Eigen::VectorXd simplexSolution;
    const int             totalIterations;
    const double          targetFunctionValue;
    const ProblemType     problemType;

public:
    SimplexResult(const Eigen::VectorXd& simplexSolution, 
                int totalIterations, 
                double targetFunctionValue, 
                ProblemType problemType)
    : simplexSolution(simplexSolution), 
        targetFunctionValue(targetFunctionValue), 
        problemType(problemType),
        totalIterations(totalIterations) {}

    ~SimplexResult() = default;

    inline const Eigen::VectorXd& getSimplexSolution() const noexcept {
        return simplexSolution;
    }

    inline int getTotalIterations() const noexcept {
        return totalIterations;
    }

    inline double getTargetFunctionValue() const noexcept {
        return targetFunctionValue;
    }

    inline ProblemType getProblemType() const noexcept {
        return problemType;
    }
};