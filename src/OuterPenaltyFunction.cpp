#include "OuterPenaltyFunction.h"

OuterPenaltyFunction::OuterPenaltyFunction(function_nd target_function) : PenaltyFunction(target_function) {
    this->equalities = std::vector<function_nd>();
}

double OuterPenaltyFunction::target_function_with_penalty(const Eigen::VectorXd& start) const {
    return target_function(start) + outer_penalty(start);
}

double OuterPenaltyFunction::outer_penalty(const Eigen::VectorXd& start) const {
    return equality_penalty(start) + inequality_penalty(start);
}

inline double OuterPenaltyFunction::equality_value_transform(double value) const {
    return std::pow(value, 2.0);
}

double OuterPenaltyFunction::equality_penalty(const Eigen::VectorXd& start) const {
    double accum = 0.0;

    for (function_nd function : equalities) {
        accum += equality_value_transform(function(start));
    }

    return accum;
}

inline double OuterPenaltyFunction::inequality_value_transform(double value) const {
    return std::pow(std::max(value, 0.0), 2.0);
}

double OuterPenaltyFunction::inequality_penalty(const Eigen::VectorXd& start) const {
    double accum = 0.0;

    for (function_nd function : inequalities) {
        accum += inequality_value_transform(function(start));
    }

    return accum;
}

size_t OuterPenaltyFunction::amount_of_equalities() const {
    return equalities.size();
}

void OuterPenaltyFunction::add_equality(function_nd equality) {
    equalities.push_back(equality);
}

void OuterPenaltyFunction::delete_equality(int index) {
    if (index < 0 || index >= equalities.size()) {
        throw std::out_of_range("Given index is out of range.");
    }

    equalities.erase(equalities.begin() + index);
}
