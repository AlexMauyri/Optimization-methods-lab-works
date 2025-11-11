#include <vector>

#include "PenaltyFunction.h"
#include "multi_dim.h"

PenaltyFunction::PenaltyFunction(function_nd target_function) {
    this->target_function = target_function;
    this->inequalities = std::vector<function_nd>();
}

size_t PenaltyFunction::amount_of_inequalities() const {
    return inequalities.size();
}

search_result_nd PenaltyFunction::compute_minimum(const Eigen::VectorXd& start) const {
    auto func = [this](const Eigen::VectorXd& x) {
        return this->target_function_with_penalty(x);
    };
    return newtone_raphson(func, start);
}

void PenaltyFunction::add_inequality(function_nd inequality) {
    inequalities.push_back(inequality);
}

void PenaltyFunction::delete_inequality(int index) {
    if (index < 0 || index >= inequalities.size()) {
        throw std::out_of_range("Given index is out of range.");
    }

    inequalities.erase(inequalities.begin() + index);
}