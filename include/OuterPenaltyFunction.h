#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>

#include "PenaltyFunction.h"
#include "multi_dim.h"
#include "search_result_nd.h"

using function_nd = std::function<double(const Eigen::VectorXd&)>;

class OuterPenaltyFunction : public PenaltyFunction {
private:
    std::vector<function_nd> equalities;

    double target_function_with_penalty(const Eigen::VectorXd& start) const override;

    double outer_penalty(const Eigen::VectorXd& start) const;

    inline double equality_value_transform(double value) const;

    double equality_penalty(const Eigen::VectorXd& start) const;

    inline double inequality_value_transform(double value) const;

    double inequality_penalty(const Eigen::VectorXd& start) const;
public:
    OuterPenaltyFunction(function_nd);

    size_t amount_of_equalities() const;

    void add_equality(function_nd equality);

    void delete_equality(int index);
};