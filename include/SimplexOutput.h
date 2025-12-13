#pragma once

#include <Eigen/Dense>
#include <memory>
#include <iomanip>
#include <type_traits>

#include "SimplexInput.h"
#include "SimplexResult.h"
#include "common.h"

class SimplexOutput final {
private:
    const bool                          isFeasible;
    const SimplexResult                 result;
    const std::shared_ptr<SimplexInput> input;

    template<typename StreamT>
    void addOriginalProblem(StreamT & output) const {
        
        output << std::string(SEPARATOR_LENGTH, '=') << "\n";
        output << "Linear programming problem\n";
        output << std::string(SEPARATOR_LENGTH, '=') << "\n\n";

        output << "Target function: ";
        if (result.getProblemType() == ProblemType::MAX) {
            output << "max F = ";
        } else {
            output << "min F = ";
        }
        addTargetFunction(output);
        
        output << "\n\nConstraints:\n\n";
        addConstraints(output);
        
        output << "\nNon-negativity conditions: ";
        for (int i = 0; i < input.get()->pricesVector.rows(); ++i) {
            output << "x" << i + 1;
            if (i < input.get()->pricesVector.rows() - 1) output << ", ";
        }
        output << " >= 0\n\n";
        
        output << std::string(SEPARATOR_LENGTH, '-') << "\n\n";
    }

    template<typename StreamT>
    void addTargetFunction(StreamT & output) const {
        bool isFirstCoefficient = true;
        for (int i = 0; i < input.get()->pricesVector.rows(); ++i) {
            double priceCoefficient = input.get()->pricesVector[i];

            if (abs(priceCoefficient) <= DOUBLE_EQUAL_PRECISION) continue;
            
            if (!isFirstCoefficient) {
                output << (priceCoefficient > 0 ? " + " : " - ");
            } else if (priceCoefficient < 0) {
                output << "-";
            }
            
            if (std::abs(std::abs(priceCoefficient) - 1.0) > DOUBLE_EQUAL_PRECISION) {
                output << std::abs(priceCoefficient);
            }
            
            output << "x" << i + 1;
            isFirstCoefficient = false;
        }
    }

    template<typename StreamT>
    void addConstraints(StreamT& output) const {
        for (int i = 0; i < input.get()->boundsVector.rows(); ++i) {
            output << "  ";
            bool isFirstCoefficient = true;
            
            for (int j = 0; j < input.get()->boundsMatrix.cols(); ++j) {
                double constraintCoefficient = input.get()->boundsMatrix(i, j);

                if (abs(constraintCoefficient) <= DOUBLE_EQUAL_PRECISION) continue;
                
                if (!isFirstCoefficient) {
                    output << (constraintCoefficient > 0 ? " + " : " - ");
                } else if (constraintCoefficient < 0) {
                    output << "-";
                }
                
                if (std::abs(std::abs(constraintCoefficient) - 1.0) > DOUBLE_EQUAL_PRECISION) {
                    output << std::abs(constraintCoefficient);
                }
                
                output << "x" << j + 1;
                isFirstCoefficient = false;
            }
            
            output << " ";
            switch (input.get()->compareSigns[i]) {
                case CompareSign::LESS_EQUAL: output << "<="; break;
                case CompareSign::EQUAL: output << "="; break;
                case CompareSign::GREATER_EQUAL: output << ">="; break;
            }
            output << " " << input.get()->boundsVector[i] << "\n";
        }
    }

    template<typename StreamT>
    void addSolution(StreamT& output) const {
        if (!isFeasible) {
            output << "Problem has no feasible solution\n\n";
            return;
        }
        
        output << "Solution result:\n\n";

        output << "Total iterations: " << result.getTotalIterations() << "\n\n";
        
        output << "Solution vector X* = (";
        for (int i = 0; i < result.getSimplexSolution().rows(); ++i) {
            output << std::fixed << std::setprecision(NUM_OF_DECIMAL_PLACES) << result.getSimplexSolution()[i];
            if (i < result.getSimplexSolution().rows() - 1) output << ", ";
        }
        output << ")\n\n";
        
        output << "Target function value at " << (result.getProblemType() == ProblemType::MAX ? "maximum" : "minimum");
        output << ": F";
        if (result.getProblemType() == ProblemType::MAX) output << "max";
        else output << "min";
        output << " = " << std::fixed << std::setprecision(NUM_OF_DECIMAL_PLACES) << result.getTargetFunctionValue() << "\n\n";
        output << std::string(SEPARATOR_LENGTH, '-') << "\n\n";
    }

    template<typename StreamT>
    void addConstraintsVerification(StreamT& output) const {  
        output << "Constraints verification:\n\n";
        for (int i = 0; i < input.get()->boundsVector.rows(); ++i) {
            double constraintValue = 0.0;
            for (int j = 0; j < result.getSimplexSolution().rows(); ++j) {
                constraintValue += input.get()->boundsMatrix(i, j) * result.getSimplexSolution()[j];
            }
            
            output << "  Constraint " << i + 1 << ": " << std::fixed << std::setprecision(NUM_OF_DECIMAL_PLACES) 
                       << constraintValue;
            
            switch (input.get()->compareSigns[i]) {
                case CompareSign::LESS_EQUAL: 
                    output << " <= " << input.get()->boundsVector[i];
                    if (constraintValue - input.get()->boundsVector[i] > DOUBLE_EQUAL_PRECISION) output << " - Violated!";
                    else output << " - Satisfied";
                    break;
                case CompareSign::EQUAL:  
                    output << " = " << input.get()->boundsVector[i];
                    if (std::abs(constraintValue - input.get()->boundsVector[i]) > DOUBLE_EQUAL_PRECISION) output << " - Violated!";
                    else output << " - Satisfied";
                    break;
                case CompareSign::GREATER_EQUAL:  
                    output << " >= " << input.get()->boundsVector[i];
                    if (input.get()->boundsVector[i] - constraintValue > DOUBLE_EQUAL_PRECISION) output << " - Violated!";
                    else output << " - Satisfied";
                    break;
            }
            output << "\n";
        }
        output << '\n' << std::string(SEPARATOR_LENGTH, '-') << "\n\n";
    }

    template<typename StreamT>
    StreamT & writeResult(StreamT& stream) const {
        addOriginalProblem(stream);
        addSolution       (stream);
        addConstraintsVerification(stream);
        return stream;    
    }



public:
    SimplexOutput(const SimplexResult& simplexResult, 
                std::shared_ptr<SimplexInput> simplexInput, 
                bool isFeasible)
    : result(simplexResult), 
    input(simplexInput), 
    isFeasible(isFeasible) {}

    SimplexOutput(const SimplexOutput& simplexOutput)
    : result(simplexOutput.result), 
    input(simplexOutput.input), 
    isFeasible(simplexOutput.isFeasible) {}

    ~SimplexOutput() = default;

    
    inline bool getIsFeasible() const noexcept {
        return isFeasible;
    }

    inline const SimplexResult& getSimplexResult() const {
        if (!isFeasible) {
            throw std::logic_error("Simplex method has no solution! Use getIsFeasible() method for checking existing.");
        }

        return result;
    }

    template<typename StreamT>
    friend StreamT& operator<<(StreamT& out, const SimplexOutput& result) {
        return result.writeResult(out);
    }
};