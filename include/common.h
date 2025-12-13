#pragma once

#include <Eigen/Dense>

#define PHI                     1.61803398874989484820
#define PSI                     0.61803398874989484820
#define N_DIM_ACCURACY          1e-5
#define N_DIM_ITERS_MAX         100
#define DX                      1e-6
#define PER_COORD_DESCEND_STEP  2.0
#define ACCURACY                1e-6
#define ITERS_MAX               50
#define DOUBLE_EQUAL_PRECISION  1e-10
#define SEPARATOR_LENGTH        60
#define NUM_OF_DECIMAL_PLACES   4

using function_nd = std::function<double(const Eigen::VectorXd&)>;

enum ProblemType {MAX, MIN};
enum CompareSign {LESS_EQUAL, EQUAL, GREATER_EQUAL};