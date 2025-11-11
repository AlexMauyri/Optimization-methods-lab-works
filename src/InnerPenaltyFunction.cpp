#include "InnerPenaltyFunction.h"

InnerPenaltyFunction::InnerPenaltyFunction(function_nd target_function) : PenaltyFunction(target_function) {}

double InnerPenaltyFunction::target_function_with_penalty(const Eigen::VectorXd& start) const {
    return target_function(start) + inner_penalty(start);
}

double InnerPenaltyFunction::inner_penalty(const Eigen::VectorXd& start) const {
    double accum = 0.0;
    
    for (function_nd function : inequalities) {
        accum += 1 / function(start);
    }

    return accum;
}