#include <Eigen/Dense>

#include "common.h"
#include "search_result_nd.h"
#include "numerics.h"

search_result_nd bisect(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional bisect method with parameters:\nleft = " << left << ";\nright = " << right 
        << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif

    search_result_nd statistic;
    statistic.type = search_method_type_nd::ND_BISECT;
    
    Eigen::VectorXd dir(left.size()), lhs(left), rhs(right);
    
    dir = direction(lhs, rhs) * 0.1 * eps;

    while (statistic.iterations != max_iterations && (statistic.accuracy = distance(lhs, rhs)) > 2.0 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic.iterations + 1 << ": left = " << lhs << "; right = " << rhs
            << "; accuracy = " << statistic.accuracy << '\n';
        #endif

        statistic.result = (lhs + rhs) * 0.5;

        Eigen::VectorXd x_l = statistic.result - dir;
        Eigen::VectorXd x_r = statistic.result + dir;

        if (function_nd(x_l) > function_nd(x_r)) {
            lhs = x_l;
        } else {
            rhs = x_r;
        }

        ++statistic.iterations;
    }

    statistic.function_probes = 2 * statistic.iterations;
    statistic.accuracy *= 0.5;

    return statistic;
}


search_result_nd golden_ratio(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional golden_ratio method with parameters:\nleft = " << left << ";\nright = " << right 
        << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif

    search_result_nd statistic;
    statistic.type = search_method_type_nd::ND_GOLDEN_RATIO;

    Eigen::VectorXd lhs(left), rhs(right);

    Eigen::VectorXd x_r = lhs + PSI * (rhs - lhs);
    Eigen::VectorXd x_l = rhs - PSI * (rhs - lhs);

    double y_r = function_nd(x_r);
    double y_l = function_nd(x_l);

    while (statistic.iterations != max_iterations && (statistic.accuracy = distance(lhs, rhs)) > 2.0 * eps) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic.iterations + 1 << ": left = " << lhs << "; right = " << rhs
            << "; accuracy = " << statistic.accuracy << '\n';
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

        ++statistic.iterations;
    }


    statistic.function_probes = statistic.iterations + 2;
    statistic.result = (lhs + rhs) * 0.5;
    statistic.accuracy *= 0.5;

    return statistic;
}


search_result_nd fibonacchi(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& left, 
    const Eigen::VectorXd& right, 
    const double eps
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional fibonacchi with parameters:\nleft = " << left << ";\nright = " << right 
        << ";\neps =" << eps << '\n';
    #endif

    search_result_nd statistic;
    statistic.type = search_method_type_nd::ND_FIBONACCHI;

    Eigen::VectorXd lhs(left), rhs(right);

    double fib_1 = 1, fib_2 = 1;
    double threshold = distance(rhs, lhs) / eps;

    while (fib_2 < threshold) {
        fib_next(fib_1, fib_2);
        ++statistic.iterations;
    }

    Eigen::VectorXd x_r = left + fib_1 * (rhs - lhs) / fib_2;
    Eigen::VectorXd x_l = left + (fib_2 - fib_1) * (rhs - lhs) / fib_2;

    double y_r = function_nd(x_r);
    double y_l = function_nd(x_l);

    for (uint64_t iterations = statistic.iterations; iterations > 0; --iterations) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic.iterations - iterations + 1 << ": left = " << lhs << "; right = " << rhs
            << "; accuracy = " << statistic.accuracy << '\n';
        #endif

        fib_prev(fib_1, fib_2);

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

    statistic.result = (lhs + rhs) * 0.5;
    statistic.accuracy = distance(rhs, lhs) * 0.5;
    statistic.function_probes = statistic.iterations + 2;

    return statistic;
}

search_result_nd per_coord_descend(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& start, 
    const double step,
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional per_coord_descend method with parameters:\nstart = " << start 
        << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif

    search_result_nd statistic;
    statistic.type = search_method_type_nd::PER_COORD_DESCEND;

    Eigen::VectorXd x_0(start), ort = Eigen::VectorXd::Zero(start.size());

    uint64_t coord_i, iteration, optimized_coord_count = 0;
    double x_i;
    for (iteration = 0; iteration < max_iterations; ++iteration) {
        #ifdef __DEBUG__
            std::cout << "Iteration #" << iteration + 1 << ": left = " << x_0 << "; right = " << x_1 << '\n';
        #endif
        coord_i = iteration % start.size();
        ort[coord_i] = 1;

        if (function_nd(x_0 - eps * ort) > function_nd(x_0 + eps * ort)) {
            ort[coord_i] = step;
        } else {
            ort[coord_i] = -step;
        }
        
        x_i = x_0[coord_i];

        auto sub_statistic = fibonacchi(function_nd, x_0, x_0 + ort, eps);

        x_0 = sub_statistic.result;
        statistic.iterations += sub_statistic.iterations;
        statistic.function_probes += sub_statistic.function_probes + 2;
    
        if (std::abs(x_0[coord_i] - x_i) < 2.0 * eps) {
            ++optimized_coord_count;
            
            if (optimized_coord_count == start.size()) {
                statistic.result = x_0;
                statistic.accuracy = sub_statistic.accuracy;
                break;
            }
        } else {
            optimized_coord_count = 0;
        }

        ort[coord_i] = 0;
    }
    
    statistic.iterations += iteration;

    return statistic;
}

search_result_nd gradient_descend(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& start, 
    const double eps, 
    const uint64_t max_iterations
) {

    #ifdef __DEBUG__
        std::cout << "Called multi dimensional gradient_descend method with parameters:\nstart = " << start 
        << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif
    search_result_nd statistic;
    statistic.type = search_method_type_nd::GRADIENT_DESCEND;

    Eigen::VectorXd curr(start.size()), prev(start), grad(start.size());

    for (; statistic.iterations != max_iterations; ++statistic.iterations) {
        grad = gradient(function_nd, prev);
        curr = prev - grad;

        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic.iterations + 1 << ": curr = " << curr << ", prev = " << prev << ", gradient = " << grad << '\n';
        #endif

        auto sub_statistic = fibonacchi(function_nd, prev, curr, eps);

        curr = sub_statistic.result;
        statistic.iterations += sub_statistic.iterations;
        statistic.function_probes += sub_statistic.function_probes + 2;

        if ((statistic.accuracy = distance(prev, curr)) < 2.0 * eps) {
            break;
        }

        prev = curr;
    }

    statistic.result = (prev + curr) * 0.5;
    statistic.accuracy *= 0.5;

    return statistic;
}

search_result_nd conj_gradient_descend(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& start, 
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional conjugate_gradient_descend method with parameters:\nstart = " << start 
        << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif
    search_result_nd statistic;
    statistic.type = search_method_type_nd::CONJ_GRADIENT_DESCEND;

    Eigen::VectorXd curr(start.size()), prev(start);
    Eigen::VectorXd curr_s, prev_s = -1.0 * gradient(function_nd, prev);

    double omega;

    for (; statistic.iterations != max_iterations; ++statistic.iterations) {
        curr = prev + prev_s;

        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic.iterations + 1 << ": curr = " << curr << ", prev = " << prev << '\n';
        #endif

        if ((statistic.accuracy = distance(prev, curr)) < 2.0 * eps) {
            break;
        }

        auto sub_statistic = fibonacchi(function_nd, prev, curr, eps);

        curr = sub_statistic.result;
        statistic.iterations += sub_statistic.iterations;
        statistic.function_probes += sub_statistic.function_probes + 2;

        curr_s = gradient(function_nd, curr);
        omega = curr_s.norm() / prev_s.norm();
        prev_s = omega * prev_s - curr_s;

        prev = curr;
    }

    statistic.result = (prev + curr) * 0.5;
    statistic.accuracy *= 0.5;

    return statistic;
}

search_result_nd newtone_raphson(
    const std::function<double(const Eigen::VectorXd)> function_nd, 
    const Eigen::VectorXd& start, 
    const double eps, 
    const uint64_t max_iterations
) {
    #ifdef __DEBUG__
        std::cout << "Called multi dimensional newtone_raphson method with parameters:\nstart = " << start 
        << ";\neps =" << eps << ";\nmax_iterations = " << max_iterations << '\n';
    #endif

    search_result_nd statistic;
    statistic.type = search_method_type_nd::NEWTONE_RAPHSON;

    Eigen::VectorXd curr, prev(start), grad;
    Eigen::MatrixXd hess(start.size(), start.size());

    for (; statistic.iterations != max_iterations; ++statistic.iterations) {
        grad = gradient(function_nd, prev);
        hess = hessian(function_nd, prev).inverse();
        curr = prev - (hess * grad);
        #ifdef __DEBUG__
            std::cout << "Iteration #" << statistic.iterations + 1 << ": curr = " << curr << ", prev = " << prev << ", gradient = " << grad << '\n';
        #endif

        if ((statistic.accuracy = distance(prev, curr)) < 2.0 * eps) {
            break;
        }

        prev = curr;
    }

    statistic.result = (prev + curr) * 0.5;
    statistic.accuracy *= 0.5;
    statistic.function_probes = statistic.iterations * 2 * start.size() * (start.size() + 2);

    return statistic;
}