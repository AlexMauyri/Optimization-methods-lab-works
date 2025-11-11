#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>

#include "PenaltyFunction.h"
#include "multi_dim.h"
#include "search_result_nd.h"

using function_nd = std::function<double(const Eigen::VectorXd&)>;

class InnerPenaltyFunction : public PenaltyFunction {
private:
    double target_function_with_penalty(const Eigen::VectorXd& start) const override;
    
    double inner_penalty(const Eigen::VectorXd& start) const;
public:
    InnerPenaltyFunction(function_nd);
};