#pragma once
#include <Eigen/Dense>
#include "common.h"
#include "search_result_nd.h"

using target_function_nd = std::function<double(const Eigen::VectorXd&)>;

search_result_nd bisect(
    target_function_nd function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps=N_DIM_ACCURACY, 
    const uint64_t max_iterations=N_DIM_ITERS_MAX
);


search_result_nd golden_ratio(
    target_function_nd function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps=N_DIM_ACCURACY, 
    const uint64_t max_iterations=N_DIM_ITERS_MAX
);


search_result_nd fibonacchi(
    target_function_nd function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps=N_DIM_ACCURACY
);

search_result_nd per_coord_descend(
    target_function_nd function_nd, 
    const Eigen::VectorXd& start,
    const double step=PER_COORD_DESCEND_STEP, 
    const double eps=N_DIM_ACCURACY, 
    const uint64_t max_iterations=N_DIM_ITERS_MAX
);

search_result_nd gradient_descend(
    target_function_nd function_nd, 
    const Eigen::VectorXd& start, 
    const double eps=N_DIM_ACCURACY, 
    const uint64_t max_iterations=N_DIM_ITERS_MAX
);

search_result_nd conj_gradient_descend(
    target_function_nd function_nd, 
    const Eigen::VectorXd& start, 
    const double eps=N_DIM_ACCURACY, 
    const uint64_t max_iterations=N_DIM_ITERS_MAX
);

search_result_nd newtone_raphson(
    target_function_nd function_nd, 
    const Eigen::VectorXd& start, 
    const double eps=N_DIM_ACCURACY, 
    const uint64_t max_iterations=N_DIM_ITERS_MAX
);