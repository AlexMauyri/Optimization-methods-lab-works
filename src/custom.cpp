#include <iostream>
#include <Eigen/Dense>

Eigen::VectorXd direction(const Eigen::VectorXd& left, const Eigen::VectorXd& right) {
    if (left.size() != right.size()) {
        throw std::runtime_error("Dimensions of vectors are not equal");
    }

    return (right - left).normalized();
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