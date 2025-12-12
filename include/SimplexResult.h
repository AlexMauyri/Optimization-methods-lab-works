#include <Eigen/Dense>
#include <iomanip>

#include "common.h"

#define SEPARATOR_LENGTH 60
#define PRECISION        6

class SimplexResult final {
private:
    const Eigen::VectorXd simplexSolution;
    const int             totalIterations;
    const double          targetFunctionValue;
    const ProblemType     problemType;

    std::stringstream     outputStream;
    std::string           formattedOutput;

    void addOriginalProblem(const Eigen::VectorXi& compareSigns,
                           const Eigen::VectorXd& pricesVector,
                           const Eigen::VectorXd& boundsVector,
                           const Eigen::MatrixXd& boundsMatrix) {
        
        outputStream << std::string(SEPARATOR_LENGTH, '=') << "\n";
        outputStream << "Linear programming problem\n";
        outputStream << std::string(SEPARATOR_LENGTH, '=') << "\n\n";

        outputStream << "Target function: ";
        if (problemType == ProblemType::MAX) {
            outputStream << "max F = ";
        } else {
            outputStream << "min F = ";
        }
        addObjectiveFunction(pricesVector);
        
        outputStream << "\n\nConstraints:\n\n";
        addConstraints(compareSigns, boundsVector, boundsMatrix);
        
        outputStream << "\nNon-negativity conditions: ";
        for (int i = 0; i < pricesVector.rows(); ++i) {
            outputStream << "x" << i + 1;
            if (i < pricesVector.rows() - 1) outputStream << ", ";
        }
        outputStream << " >= 0\n\n";
        
        outputStream << std::string(SEPARATOR_LENGTH, '-') << "\n\n";
    }

    void addObjectiveFunction(const Eigen::VectorXd& pricesVector) {
        bool isFirstCoefficient = true;
        for (int i = 0; i < pricesVector.rows(); ++i) {
            double priceCoefficient = pricesVector[i];

            if (abs(priceCoefficient) <= DOUBLE_EQUAL_PRECISION) continue;
            
            if (!isFirstCoefficient) {
                outputStream << (priceCoefficient > 0 ? " + " : " - ");
            } else if (priceCoefficient < 0) {
                outputStream << "-";
            }
            
            if (std::abs(std::abs(priceCoefficient) - 1.0) > DOUBLE_EQUAL_PRECISION) {
                outputStream << std::abs(priceCoefficient);
            }
            
            outputStream << "x" << i + 1;
            isFirstCoefficient = false;
        }
    }

    void addConstraints(const Eigen::VectorXi& compareSigns,
                       const Eigen::VectorXd& boundsVector,
                       const Eigen::MatrixXd& boundsMatrix) {
        for (int i = 0; i < boundsVector.rows(); ++i) {
            outputStream << "  ";
            bool isFirstCoefficient = true;
            
            for (int j = 0; j < boundsMatrix.cols(); ++j) {
                double constraintCoefficient = boundsMatrix(i, j);

                if (abs(constraintCoefficient) <= DOUBLE_EQUAL_PRECISION) continue;
                
                if (!isFirstCoefficient) {
                    outputStream << (constraintCoefficient > 0 ? " + " : " - ");
                } else if (constraintCoefficient < 0) {
                    outputStream << "-";
                }
                
                if (std::abs(std::abs(constraintCoefficient) - 1.0) > DOUBLE_EQUAL_PRECISION) {
                    outputStream << std::abs(constraintCoefficient);
                }
                
                outputStream << "x" << j + 1;
                isFirstCoefficient = false;
            }
            
            outputStream << " ";
            switch (compareSigns[i]) {
                case CompareSign::LESS_EQUAL: outputStream << "<="; break;
                case CompareSign::EQUAL: outputStream << "="; break;
                case CompareSign::GREATER_EQUAL: outputStream << ">="; break;
            }
            outputStream << " " << boundsVector[i] << "\n";
        }
    }

    void addSolution() {
        if (simplexSolution.size() == 0) {
            outputStream << "Problem has no feasible solution\n";
            return;
        }
        
        outputStream << "Solution result:\n\n";

        outputStream << "Total iterations: " << totalIterations << "\n\n";
        
        outputStream << "Solution vector X* = (";
        for (int i = 0; i < simplexSolution.rows(); ++i) {
            outputStream << std::fixed << std::setprecision(PRECISION) << simplexSolution[i];
            if (i < simplexSolution.rows() - 1) outputStream << ", ";
        }
        outputStream << ")\n\n";
        
        outputStream << "Target function value at " << (problemType == ProblemType::MAX ? "maximum" : "minimum");
        outputStream << ": F";
        if (problemType == ProblemType::MAX) outputStream << "max";
        else outputStream << "min";
        outputStream << " = " << std::fixed << std::setprecision(PRECISION) << targetFunctionValue << "\n\n";
        outputStream << std::string(SEPARATOR_LENGTH, '-') << "\n\n";
    }

    void addConstraintsVerification(const Eigen::VectorXi& compareSigns,
                                   const Eigen::VectorXd& boundsVector,
                                   const Eigen::MatrixXd& boundsMatrix) {  
        outputStream << "Constraints verification:\n\n";
        for (int i = 0; i < boundsVector.rows(); ++i) {
            double constraintValue = 0.0;
            for (int j = 0; j < simplexSolution.rows(); ++j) {
                constraintValue += boundsMatrix(i, j) * simplexSolution[j];
            }
            
            outputStream << "  Constraint " << i + 1 << ": " << std::fixed << std::setprecision(PRECISION) 
                       << constraintValue;
            
            switch (compareSigns[i]) {
                case CompareSign::LESS_EQUAL: 
                    outputStream << " <= " << boundsVector[i];
                    if (constraintValue - boundsVector[i] > 1e-6) outputStream << " - Violated!";
                    else outputStream << " - Satisfied";
                    break;
                case CompareSign::EQUAL:  
                    outputStream << " = " << boundsVector[i];
                    if (std::abs(constraintValue - boundsVector[i]) > 1e-6) outputStream << " - Violated!";
                    else outputStream << " - Satisfied";
                    break;
                case CompareSign::GREATER_EQUAL:  
                    outputStream << " >= " << boundsVector[i];
                    if (boundsVector[i] - constraintValue > 1e-6) outputStream << " - Violated!";
                    else outputStream << " - Satisfied";
                    break;
            }
            outputStream << "\n";
        }
        outputStream << '\n' << std::string(SEPARATOR_LENGTH, '-') << "\n\n";
    }

    void finalizeOutput() {
        formattedOutput = outputStream.str();
        outputStream.clear();
        outputStream.str(std::string());
    }

public:
    SimplexResult(const Eigen::VectorXd& simplexSolution, double targetFunctionValue, ProblemType problemType, int totalIterations)
        : simplexSolution(simplexSolution), 
        targetFunctionValue(targetFunctionValue), 
        problemType(problemType),
        totalIterations(totalIterations) {}

    void writeResult(const Eigen::VectorXi& compareSigns,
                    const Eigen::VectorXd& pricesVector,
                    const Eigen::VectorXd& boundsVector,
                    const Eigen::MatrixXd& boundsMatrix) {
        addOriginalProblem(compareSigns, pricesVector, boundsVector, boundsMatrix);
        addSolution();
        addConstraintsVerification(compareSigns, boundsVector, boundsMatrix);
        finalizeOutput();
    }

    friend std::ostream& operator<<(std::ostream& os, const SimplexResult& result) {
        os << result.formattedOutput;
        return os;
    }
};