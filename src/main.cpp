#include <functional>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <Eigen/Dense>

#include "one_dim.h"
#include "multi_dim.h"
#include "common.h"
#include "numerics.h"
#include "InnerPenaltyFunction.h"
#include "OuterPenaltyFunction.h"
#include "Simplex.h"

double test_func(const double x) {
	return (x - 1) * (x - 5); //3
}

double test_func_2(const Eigen::VectorXd& x) {
	return (x[0] - 5) * x[0] + (x[1] - 3) * x[1]; //[2.5, 1.5]
}

double psi1(const Eigen::VectorXd& x) {
	return 5 - x[0] * 2 + x[1] * 3;
}

double psi2(const Eigen::VectorXd& x) {
	return 6 + x[0] * 3 - x[1];
}

double phi1(const Eigen::VectorXd& x) {
	return x[0] - x[1] - 1;
}

void lab1(const std::function<double(double)> function) {
    double x0 = 0;
    double x1 = 10;

    std::cout << bisect(function, x0, x1) << '\n';
    std::cout << golden_ratio(function, x0, x1) << '\n';
    std::cout << fibonacchi(function, x0, x1, ACCURACY + 2e-7) << '\n';
}

void lab2(function_nd function_nd) {
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

void lab3(function_nd function_nd) {
    Eigen::VectorXd start(2);
    start << -14, -33.98;

    std::cout << gradient_descend(function_nd, start) << '\n';
    std::cout << conj_gradient_descend(function_nd, start) << '\n';
    std::cout << newton(function_nd, start) << '\n';

    start[0] = -12.0;
    start[1] = -15.0;
   
    std::cout << newton(function_nd, start) << '\n';

    InnerPenaltyFunction penalty_function1 = InnerPenaltyFunction(test_func_2);
    penalty_function1.addInequality(psi1);
    penalty_function1.addInequality(psi2);

    std::cout << penalty_function1.computeMinimum(start) << '\n';

    OuterPenaltyFunction penalty_function2 = OuterPenaltyFunction(test_func_2);
    penalty_function2.addInequality(psi1);
    penalty_function2.addInequality(psi2);
    penalty_function2.addEquality(phi1);

    std::cout << penalty_function2.computeMinimum(start) << '\n';
}

void lab4() {
    Eigen::VectorXi compareSigns(3);
    Eigen::VectorXd pricesVector(2);
    Eigen::VectorXd boundsVector(3);
    Eigen::MatrixXd boundsMatrix(3, 2);

    compareSigns << 2, 2, 2;
    pricesVector << 2, 3;
    boundsVector << 40, 28, 14;
    boundsMatrix << -2, 6, 3, 2, 2, -1;

    std::cout << "Starting Simplex method\n";

    Simplex simplex {boundsMatrix, boundsVector, pricesVector, compareSigns};

    std::cout << simplex.solve(ProblemType::MIN);

    Eigen::VectorXi compareSigns1(2);
    Eigen::VectorXd pricesVector1(3);
    Eigen::VectorXd boundsVector1(2);
    Eigen::MatrixXd boundsMatrix1(2, 3);

    compareSigns1 << 2, 0;
    pricesVector1 << 1, 0, 0;
    boundsVector1 << 2, 20;
    boundsMatrix1 << 1, 1, 1, 4, 5, 5;

    Simplex simplex1 {boundsMatrix1, boundsVector1, pricesVector1, compareSigns1};

    std::cout << simplex1.solve(ProblemType::MAX);
}

int main() {
    //lab1(test_func);
    //lab2(test_func_2);
    //lab3(test_func_2);
    lab4();

    return 0;
}