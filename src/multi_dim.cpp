#include <Eigen/Dense>
#include <tuple>
#include "common.h"
#include "search_result_nd.h"
#include "custom.h"

search_result_nd* bisect(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional bisect method with parameters:\nleft = ";
        custom_vector_print(std::cout, left); 
        std::cout << ";\nright = ";
        custom_vector_print(std::cout, right); 
        std::cout << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif

    search_result_nd* statistic = new search_result_nd();
    statistic->type = search_method_type_nd::ND_BISECT;
    
    Eigen::VectorXd dir(left.size()), lhs(left), rhs(right);
    
    dir = direction(lhs, rhs) * 0.1 * eps;

    while (statistic->iterations != max_iterations && (statistic->accuracy = distance(lhs, rhs)) > 2 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations + 1 << ": left = ";
            custom_vector_print(std::cout, lhs);  
            std::cout << "; right = "; 
            custom_vector_print(std::cout, rhs);
            std::cout << '\n';
        #endif

        statistic->result = (lhs + rhs) * 0.5;

        Eigen::VectorXd x_l = statistic->result - dir;
        Eigen::VectorXd x_r = statistic->result + dir;

        if (function_nd(x_l) > function_nd(x_r)) {
            lhs = x_l;
        } else {
            rhs = x_r;
        }

        ++statistic->iterations;
    }

    statistic->function_probes = 2 * statistic->iterations;

    return statistic;
}


search_result_nd* golden_ratio(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional golden_ratio method with parameters:\nleft = ";
        custom_vector_print(std::cout, left); 
        std::cout << ";\nright = ";
        custom_vector_print(std::cout, right); 
        std::cout << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif

    search_result_nd* statistic = new search_result_nd();
    statistic->type = search_method_type_nd::ND_GOLDEN_RATIO;

    Eigen::VectorXd lhs(left), rhs(right);

    Eigen::VectorXd x_r = lhs + PSI * (rhs - lhs);
    Eigen::VectorXd x_l = rhs - PSI * (rhs - lhs);

    double y_r = function_nd(x_r);
    double y_l = function_nd(x_l);

    while (statistic->iterations != max_iterations && (statistic->accuracy = distance(lhs, rhs)) > 2 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations + 1 << ": left = ";
            custom_vector_print(std::cout, lhs);  
            std::cout << "; right = "; 
            custom_vector_print(std::cout, rhs);
            std::cout << '\n';
        #endif

        if (y_l > y_r) {
            lhs = x_l;
            x_l = x_r;
            y_l = y_r;
            x_r = lhs + PSI * (rhs - lhs);
            y_r = function_nd(x_r);
        } else {
            rhs = x_r;
            x_r = x_l;
            y_r = y_l;
            x_l = rhs - PSI * (rhs - lhs);
            y_l = function_nd(x_l);
        }

        ++statistic->iterations;
    }


    statistic->function_probes = statistic->iterations + 2;
    statistic->result = (lhs + rhs) * 0.5;

    return statistic;
}


search_result_nd* fibonacchi(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional fibonacchi with parameters:\nleft = ";
        custom_vector_print(std::cout, left); 
        std::cout << ";\nright = ";
        custom_vector_print(std::cout, right); 
        std::cout << ";\neps =" << eps << '\n';
    #endif

    search_result_nd* statistic = new search_result_nd();
    statistic->type = search_method_type_nd::ND_FIBONACCHI;

    Eigen::VectorXd lhs(left), rhs(right);

    uint64_t fib_temp = 0, fib_1 = 1, fib_2 = 1;
    double threshold = distance(rhs, lhs) / eps;

    while (fib_2 < threshold) {
        fib_temp = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_temp;
        ++statistic->iterations;
    }

    Eigen::VectorXd x_r = left + static_cast<double>(fib_1) / fib_2 * (rhs - lhs);
    Eigen::VectorXd x_l = left + static_cast<double>(fib_2 - fib_1) / fib_2 * (rhs - lhs);

    double y_r = function_nd(x_r);
    double y_l = function_nd(x_l);

    for (uint64_t iterations = statistic->iterations; iterations > 0; --iterations) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations - iterations + 1 << ": left = ";
            custom_vector_print(std::cout, lhs);  
            std::cout << "; right = "; 
            custom_vector_print(std::cout, rhs);
            std::cout << '\n';
        #endif

        fib_temp = fib_2 - fib_1;
        fib_2 = fib_1;
        fib_1 = fib_temp;

        if (y_l > y_r) {
            lhs = x_l;
            x_l = x_r;
            y_l = y_r;
            x_r = lhs + static_cast<double>(fib_1) / fib_2 * (rhs - lhs);
            y_r = function_nd(x_r);
        } else {
            rhs = x_r;
            x_r = x_l;
            y_r = y_l;
            x_l = lhs + static_cast<double>(fib_2 - fib_1) / fib_2 * (rhs - lhs);
            y_l = function_nd(x_l);
        }
    }

    statistic->result = (lhs + rhs) * 0.5;
    statistic->accuracy = distance(rhs, lhs) / 2;
    statistic->function_probes = statistic->iterations + 2;

    return statistic;
}

search_result_nd* per_coord_descend(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& start, 
    const double step,
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional per_coord_descend method with parameters:\nstart = ";
        custom_vector_print(std::cout, start); 
        std::cout << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif

    search_result_nd* statistic = new search_result_nd();
    statistic->type = search_method_type_nd::PER_COORD_DESCEND;

    Eigen::VectorXd x_0(start), x_1(start);

    uint64_t coord_i, iteration, optimized_coord_count = 0;
    double x_i;
    for (iteration = 0; iteration < max_iterations; ++iteration) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << iteration + 1 << ": left = ";
            custom_vector_print(std::cout, x_0);  
            std::cout << "; right = "; 
            custom_vector_print(std::cout, x_1);
            std::cout << '\n';
        #endif
        coord_i = iteration % start.size();

        x_1[coord_i] -= eps;
        double y_0 = function_nd(x_1);
        x_1[coord_i] += 2.0 * eps;
        double y_1 = function_nd(x_1);
        x_1[coord_i] -= eps;

        if (y_0 > y_1) {
            x_1[coord_i] += step;
        } else {
            x_1[coord_i] -= step;
        }
        
        x_i = x_0[coord_i];

        std::cout.setstate(std::ios_base::failbit);
        auto sub_statistic = fibonacchi(function_nd, x_0, x_1, eps);
        std::cout.clear();

        statistic->result = sub_statistic->result;
        statistic->accuracy = sub_statistic->accuracy;
        statistic->function_probes += sub_statistic->function_probes + 2;

        x_0 = sub_statistic->result;
    
        if (std::abs(x_1[coord_i] - x_i) < 2.0 * eps) {
            ++optimized_coord_count;
            
            if (optimized_coord_count == x_1.size()) {
                break;
            }
        }
    }
    statistic->iterations = iteration;

    return statistic;
}

search_result_nd* gradient_descend(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& start, 
    const double eps=N_DIM_ACCURACY, 
    const uint64_t max_iterations=N_DIM_ITERS_MAX
) {

    #ifdef __DEBUG__
        std::cout << "Called multi dimensional gradient_descend method with parameters:\nstart = ";
        custom_vector_print(std::cout, start); 
        std::cout << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif
    search_result_nd* statistic = new search_result_nd();
    statistic->type = search_method_type_nd::GRADIENT_DESCEND;

    Eigen::VectorXd x_i(start), x_i_1(start.size()), grad(start.size());

    do {
        grad = gradient(function_nd, x_i);
        x_i_1 = x_i;
        x_i -= grad;

        std::cout.setstate(std::ios_base::failbit);
        auto sub_statistic = fibonacchi(function_nd, x_i_1, x_i, eps);
        std::cout.clear();

        x_i = sub_statistic->result;
        statistic->function_probes += sub_statistic->function_probes + 2;

        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic->iterations + 1 << ": x_i = ";
            custom_vector_print(std::cout, x_i);
            std::cout << ", x_i-1 = ";
            custom_vector_print(std::cout, x_i_1);
            std::cout << ", gradient = ";
            custom_vector_print(std::cout, grad);
            std::cout << '\n';
        #endif

        ++statistic->iterations;
    } while (statistic->iterations != max_iterations && (statistic->accuracy = distance(x_i, x_i_1)) >= 2.0 * eps);

    statistic->result = (x_i_1 + x_i) * 0.5;

    return statistic;
}