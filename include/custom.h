#include <Eigen/Dense>

Eigen::VectorXd direction(const Eigen::VectorXd& left, const Eigen::VectorXd& right);

Eigen::VectorXd gradient(const std::function<double(const Eigen::VectorXd)> function_nd, const Eigen::VectorXd& point);

double distance(const Eigen::VectorXd& left, const Eigen::VectorXd& right);

void custom_vector_print(std::ostream& out, const Eigen::VectorXd& vec);

