#include <functional>
#include <iostream>
#include <Eigen/Dense>

#include "one_dim.h"
#include "multi_dim.h"
#include "common.h"
#include "numerics.h"

double test_func(const double x)
{
	return (x - 1) * (x - 5); //3
}

double test_func_2(const Eigen::VectorXd& x)
{
	return (x[0] - 5) * x[0] + (x[1] - 3) * x[1]; //[2.5, 1.5]
}

void lab1(const std::function<double(double)> function) {
    double x0 = 0;
    double x1 = 10;

    std::cout << bisect(function, x0, x1) << '\n';
    std::cout << golden_ratio(function, x0, x1) << '\n';
    std::cout << fibonacchi(function, x0, x1, ACCURACY + 2e-7) << '\n';
}

void lab2(std::function<double(const Eigen::VectorXd&)> function_nd) {
    Eigen::VectorXd x_0(2);
	Eigen::VectorXd x_1(2);

    x_0 << 0.0, 0.0;
    x_1 << 5.0, 3.0;
    
	std::cout << bisect(function_nd, x_1, x_0) << '\n';
	std::cout << golden_ratio(function_nd, x_1, x_0) << '\n';
	std::cout << fibonacchi(function_nd, x_1, x_0) << '\n';
	
    Eigen::VectorXd start(2);
    start << -14, -33.98;
	std::cout << per_coord_descend(function_nd, start) << '\n';
}

void lab3(std::function<double(const Eigen::VectorXd&)> function_nd) {
    Eigen::VectorXd start(2);
    start << -14, -33.98;

    std::cout << gradient_descend(function_nd, start) << '\n';
    std::cout << conj_gradient_descend(function_nd, start) << '\n';
    std::cout << newtone_raphson(function_nd, start) << '\n';
}

int main() {
    //lab1(test_func);
    lab2(test_func_2);
    //lab3(test_func_2);
    

    return 0;
}