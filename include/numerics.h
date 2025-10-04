#include <Eigen/Dense>

Eigen::VectorXd direction(const Eigen::VectorXd& left, const Eigen::VectorXd& right);

double partial(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point, uint32_t index);

double partial2(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point, uint32_t index1, uint32_t index2);

Eigen::VectorXd gradient(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point);

Eigen::MatrixXd hessian(std::function<double(const Eigen::VectorXd)> function_nd, Eigen::VectorXd& point);

double distance(const Eigen::VectorXd& left, const Eigen::VectorXd& right);

void custom_vector_print(std::ostream& out, const Eigen::VectorXd& vec);

std::ostream& operator<<(std::ostream& stream, const Eigen::VectorXd& vec);

inline void fib_next(double& fib_1, double& fib_2) {
    const double fib_temp(fib_1);
    fib_1 = fib_2;
    fib_2 += fib_temp;
}

inline void fib_prev(double& fib_1, double& fib_2) {
    const double fib_temp(fib_2 - fib_1);
    fib_2 = fib_1;
    fib_1 = fib_temp;
}