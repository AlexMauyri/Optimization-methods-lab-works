#include <iostream>
#include <Eigen/Dense>

#include "common.h"

Eigen::VectorXd direction(const Eigen::VectorXd& left, const Eigen::VectorXd& right) {
    if (left.size() != right.size()) {
        throw std::runtime_error("Dimensions of vectors are not equal");
    }

    return (right - left).normalized();
}

Eigen::VectorXd gradient(const std::function<double(const Eigen::VectorXd)> function_nd, const Eigen::VectorXd& point) {
    Eigen::VectorXd grad(point.size()), copy(point);

    for (int i = 0; i < grad.size(); ++i) {
        
        copy[i] += DX;
        double y_1 = function_nd(copy);
        copy[i] -= 2.0 * DX;
        double y_2 = function_nd(copy);

        grad[i] = ((y_1 - y_2) / (2.0 * DX));
    }

    return grad;
}

double distance(const Eigen::VectorXd& left, const Eigen::VectorXd& right) {
    if (left.size() != right.size()) {
        throw std::runtime_error("Dimensions of vectors are not equal");
    }

    return (right - left).norm();
}

void custom_vector_print(std::ostream& stream, const Eigen::VectorXd& vec) {
    stream << '[' << vec[0];
    for (int i = 1; i < vec.size(); ++i) {
        stream << ',' << ' ' << vec[i]; 
    }
    stream << ']';
}