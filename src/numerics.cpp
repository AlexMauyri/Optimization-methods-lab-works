#include <iostream>
#include <Eigen/Dense>

#include "common.h"
#include "numerics.h"

Eigen::VectorXd direction(const Eigen::VectorXd& left, const Eigen::VectorXd& right) {
    if (left.size() != right.size()) {
        throw std::runtime_error("Dimensions of vectors are not equal");
    }

    return (right - left).normalized();
}

double partial(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point, uint32_t index) {
    if (index >= point.size()) {
        throw std::runtime_error("Index value is out of vector indices range");
    }

    point[index] += DX;
    double y_1 = function_nd(point);
    point[index] -= 2.0 * DX;
    double y_2 = function_nd(point);
    point[index] += DX;

    return ((y_1 - y_2) / (2.0 * DX));
}

double partial2(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point, uint32_t index1, uint32_t index2) {
    if (index2 >= point.size()) {
        throw std::runtime_error("Index value is out of vector indices range");
    }
    point[index2] += DX;
    double y_1 = partial(function_nd, point, index1);
    point[index2] -= 2.0 * DX;
    double y_2 = partial(function_nd, point, index1);
    point[index2] += DX;

    return ((y_1 - y_2) / (2.0 * DX));
}

Eigen::VectorXd gradient(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point) {
    Eigen::VectorXd grad(point.size());

    for (int i = 0; i < grad.size(); ++i) {
        grad[i] = partial(function_nd, point, i);
    }

    return grad;
}

Eigen::MatrixXd hessian(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point) {
    Eigen::MatrixXd result(point.size(), point.size());

    for (uint32_t row = 0; row < result.rows(); ++row) {
        for (uint32_t col = row; col < result.cols(); ++col) {
            result(col, row) = result(row, col) = partial2(function_nd, point, row, col);
        }
    }

    return result;
}

double distance(const Eigen::VectorXd& left, const Eigen::VectorXd& right) {
    if (left.size() != right.size()) {
        throw std::runtime_error("Dimensions of vectors are not equal");
    }

    return (right - left).norm();
}

std::ostream& operator<<(std::ostream& stream, const Eigen::VectorXd& vec) {
    stream << '[' << vec[0];
    for (int i = 1; i < vec.size(); ++i) {
        stream << ',' << ' ' << vec[i]; 
    }
    stream << ']';
}