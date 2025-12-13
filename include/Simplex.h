#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <memory>

#include "SimplexInput.h"
#include "SimplexOutput.h"
#include "common.h"

enum class SimplexStepCode {NEXT_STEP, NO_MAIN_COLUMN, NO_MAIN_ROW};

class Simplex final {
private:
    static constexpr int                SIMPLEX_MAX_ITERATIONS = 50;
    static constexpr double             SIMPLEX_ACCURACY       = 1e-6;

    const std::shared_ptr<SimplexInput> input;
    Eigen::VectorXi                     basicVariableIndices;
    Eigen::VectorXi                     naturalVariableIndices;
    Eigen::VectorXi                     virtualVariableIndices;
    Eigen::MatrixXd                     simplexTable;
    ProblemType                         problemType;

    void buildSimplexTable() {
        allocateMemory();
        makeCoefficientsPositive();
        initializeAllVariables();
        populateConstraintRows();
        addTargetFunctionRow();
        if (!isTargetFunctionModified()) return;
        addArtificialTargetFunctionRow();
    }

    inline void allocateMemory() {
        size_t naturalNumber = input.get()->pricesVector.rows() + input.get()->boundsVector.rows(), basicNumber = 0, virtualNumber = 0;

        for (size_t index = 0; index < input.get()->compareSigns.rows(); ++index) {
            if ((input.get()->compareSigns[index] == CompareSign::LESS_EQUAL && input.get()->boundsVector[index] > 0) 
                || (input.get()->compareSigns[index] == CompareSign::GREATER_EQUAL && input.get()->boundsVector[index] < 0)) {
                basicNumber += 1;
            } else {
                basicNumber += 2;
                virtualNumber += 1;
            }
        }

        basicVariableIndices = Eigen::VectorXi(basicNumber);
        naturalVariableIndices = Eigen::VectorXi(naturalNumber);
        virtualVariableIndices = Eigen::VectorXi(virtualNumber);

        #ifdef __DEBUG__
            std::cout << "Allocated vectors\n";
            std::cout << "Natural variables amount: " << naturalNumber << "\n";
            std::cout << "Basic variables amount: " << basicNumber << "\n";
            std::cout << "Virtual variables amount: " << virtualNumber << "\n";
        #endif

        simplexTable = Eigen::MatrixXd::Zero(input.get()->boundsVector.rows() + 1 + isTargetFunctionModified(), naturalNumber + virtualNumber + 1);

        #ifdef __DEBUG__
            std::cout << "Allocated simplex table\n";
            std::cout << "Rows: " << simplexTable.rows() << "\n";
            std::cout << "Cols: " << simplexTable.cols() << "\n";
        #endif
    }

    inline void makeCoefficientsPositive() noexcept {
        for (size_t rowIndex = 0; rowIndex < input.get()->boundsMatrix.rows(); ++rowIndex) {
            if (input.get()->boundsVector[rowIndex] >= 0) continue;
            input.get()->boundsVector[rowIndex] *= -1.0;
            simplexTable.row(rowIndex) *= -1.0;
            
            if (input.get()->compareSigns[rowIndex] == CompareSign::EQUAL) continue;
            input.get()->compareSigns[rowIndex] = input.get()->compareSigns[rowIndex] == CompareSign::LESS_EQUAL? CompareSign::GREATER_EQUAL : CompareSign::LESS_EQUAL;
        }

        #ifdef __DEBUG__
            std::cout << "Made coefficients positive\n";
        #endif
    }

    inline void initializeAllVariables() noexcept {
        size_t naturalVariablesCount = 0;
        size_t basicVariablesCount = 0;
        size_t virtualVariablesCount = 0;

        for (size_t index = 0; index < input.get()->pricesVector.rows(); ++index) {
            naturalVariableIndices[naturalVariablesCount++] = index;
        }

        #ifdef __DEBUG__
            std::cout << "Added variables from target function\n";
        #endif

        size_t lastFilledColIndex = input.get()->boundsMatrix.cols() - 1;

        int additionalVariableIndex, virtualVariableIndex;

        for (size_t inequalityIndex = 0; inequalityIndex < input.get()->boundsMatrix.rows(); ++inequalityIndex) {
            createArtificialVariables(inequalityIndex, additionalVariableIndex, virtualVariableIndex, lastFilledColIndex);

            naturalVariableIndices[naturalVariablesCount++] = additionalVariableIndex;

            if (virtualVariableIndex != -1) {
                basicVariableIndices[basicVariablesCount++] = virtualVariableIndex;
                virtualVariableIndices[virtualVariablesCount++] = virtualVariableIndex;
            }

            basicVariableIndices[basicVariablesCount++] = additionalVariableIndex;
        }

        #ifdef __DEBUG__
            std::cout << "Created additional and virtual variables\n";
            std::cout << "Natural variable indices: " << naturalVariableIndices << '\n';
            std::cout << "Basic variable indices: " << basicVariableIndices << '\n';
            std::cout << "Virtual variable indices: " << virtualVariableIndices << '\n';
        #endif
    }

    inline void createArtificialVariables(int inequalityIndex, 
                                int& additionalVariableIndex, 
                                int& virtualVariableIndex,
                                size_t& lastFilledColIndex) noexcept {
        if (input.get()->compareSigns[inequalityIndex] == CompareSign::GREATER_EQUAL) {
            simplexTable(inequalityIndex, ++lastFilledColIndex) = -1.0;
            simplexTable(inequalityIndex, ++lastFilledColIndex) = 1.0;
            additionalVariableIndex = lastFilledColIndex - 1;
            virtualVariableIndex = lastFilledColIndex;
        } else {
            simplexTable(inequalityIndex, ++lastFilledColIndex) = 1.0;
            additionalVariableIndex = lastFilledColIndex;
            virtualVariableIndex = input.get()->compareSigns[inequalityIndex] == CompareSign::EQUAL? additionalVariableIndex : -1;
        }
    }

    inline void populateConstraintRows() noexcept {
        for (size_t rowIndex = 0; rowIndex < input.get()->boundsMatrix.rows(); ++rowIndex) {
            for (size_t colIndex = 0; colIndex < input.get()->boundsMatrix.cols(); ++colIndex) {
                simplexTable(rowIndex, colIndex) = input.get()->boundsMatrix(rowIndex, colIndex);
            }
        }

        #ifdef __DEBUG__
            std::cout << "Carried values from bounds matrix to simplex table\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif

        simplexTable.block(0, simplexTable.cols() - 1, input.get()->boundsVector.rows(), 1) = input.get()->boundsVector;

        #ifdef __DEBUG__
            std::cout << "Carried free coefficients from bounds vector to simplex table\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif
    }

    inline void addTargetFunctionRow() noexcept {
        for (size_t index = 0; index < input.get()->pricesVector.rows(); ++index) {
            double sign = problemType == ProblemType::MAX ? -1.0 : 1.0;
            simplexTable(input.get()->boundsMatrix.rows(), index) = sign * input.get()->pricesVector[index];
        }

        #ifdef __DEBUG__
            std::cout << "Adding target function row\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif
    }

    inline void addArtificialTargetFunctionRow() noexcept {
        for (size_t virtualVariableIndex : virtualVariableIndices) {
            simplexTable(input.get()->boundsMatrix.rows() + 1, virtualVariableIndex) = 1.0;
        }

        #ifdef __DEBUG__
            std::cout << "Added artificial target function row\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif
    }

    inline void excludeVirtualVariables() noexcept {
        if (!isTargetFunctionModified()) return;

        const size_t lastRowId = simplexTable.rows() - 1;

        for (size_t virtualVariableIndex : virtualVariableIndices) {
            for (size_t rowIndex = 0; rowIndex < simplexTable.rows(); ++rowIndex) {
                if (simplexTable(rowIndex, virtualVariableIndex) == 0) continue;
                simplexTable.row(lastRowId) -= simplexTable.row(rowIndex);
                break;
            }
        }

        #ifdef __DEBUG__
            std::cout << "Excluded virtual variables\n";
        #endif
    }

    int getMainColumn() const noexcept {
        Eigen::VectorXd row = simplexTable.row(simplexTable.rows() - 1);
        double minimumNegative = 0;
        int mainCol = -1;

        for (size_t index = 0; index < row.rows() - 1; ++index) {
	        if (row[index] >= minimumNegative) continue;
            minimumNegative = row[index];
            mainCol = index;
        }

        if (!isTargetFunctionModified() || mainCol != -1) return mainCol;

        row = simplexTable.row(simplexTable.rows() - 2);

        for (size_t index : naturalVariableIndices) {
            if (row[index] >= minimumNegative) continue;
            minimumNegative = row[index];
            mainCol = index;
        }

        return mainCol;
    }

    int getMainRow(int mainCol) const noexcept {
        double minimumValue = std::numeric_limits<double>::max(), mainElement;
        int mainRow = -1, lastColIndex = simplexTable.cols() - 1;
        int targetFunctionRowIndex = isTargetFunctionModified()? simplexTable.rows() - 2 : simplexTable.rows() - 1;

        for (size_t index = 0; index < targetFunctionRowIndex; ++index) {
            mainElement = simplexTable(index, mainCol);
            if (mainElement <= 0.0) continue;
            if (simplexTable(index, lastColIndex) / mainElement > minimumValue) continue;

            minimumValue = simplexTable(index, lastColIndex) / mainElement;
            mainRow = index; 
        }

        return mainRow;
    }

    bool validateSolution() const noexcept {
        double targetFunctionValue = 0.0;

        const size_t targetFunctionRowIndex = isTargetFunctionModified()? simplexTable.rows() - 2 : simplexTable.rows() - 1;
        const size_t targetFunctionColIndex = simplexTable.cols() - 1;

        for (size_t index = 0; index < basicVariableIndices.rows(); ++index) {
            if (basicVariableIndices[index] >= input.get()->pricesVector.rows()) {
                continue;
            }

            targetFunctionValue += simplexTable(index, targetFunctionColIndex) * input.get()->pricesVector[basicVariableIndices[index]];
        }

        double sign = problemType == ProblemType::MAX? -1.0 : 1.0;
        double targetFunctionTableValue = simplexTable(targetFunctionRowIndex, targetFunctionColIndex);
        #ifdef __DEBUG__
            std::cout << "Value of target function: " << targetFunctionValue << "; Value of target function in table: " << targetFunctionTableValue << '\n';
        #endif
        if (std::abs(targetFunctionValue + sign * targetFunctionTableValue) < SIMPLEX_ACCURACY) {
            if (isTargetFunctionModified()) {
                return std::abs(simplexTable(simplexTable.rows() - 1, targetFunctionColIndex)) < SIMPLEX_ACCURACY;
            }
            return true;
        }

        return false;
    }

    Eigen::VectorXd getSimplexSolution(bool onlyNaturalArgs) const noexcept {
        Eigen::VectorXd solution = Eigen::VectorXd::Zero(onlyNaturalArgs? input.get()->pricesVector.rows() : simplexTable.cols() - 1);

        for (size_t index = 0; index < basicVariableIndices.rows(); ++index) {
            if (basicVariableIndices[index] >= solution.rows()) continue;
            solution[basicVariableIndices[index]] = simplexTable(index, simplexTable.cols() - 1);
        }

        return solution;
    }

    SimplexStepCode simplexStep() noexcept {
        double mainElement;
        int mainRow, mainCol;

        mainCol = getMainColumn();

        if (mainCol == -1) {
            return SimplexStepCode::NO_MAIN_COLUMN;
        }

        mainRow = getMainRow(mainCol);

        if (mainRow == -1) {
            return SimplexStepCode::NO_MAIN_ROW;
        }

        basicVariableIndices[mainRow] = mainCol;
        mainElement = simplexTable(mainRow, mainCol);
        simplexTable.row(mainRow) /= mainElement;

        for (size_t index = 0; index < simplexTable.rows(); ++index) {
            if (index == mainRow) continue;
            simplexTable.row(index) -= simplexTable(index, mainCol) * simplexTable.row(mainRow);
        }

        return SimplexStepCode::NEXT_STEP;
    }

    SimplexOutput generateOutput(const int iteration) const noexcept {
        const bool isFeasible = validateSolution();
        const Eigen::VectorXd& simplexSolution = isFeasible? getSimplexSolution(true) : Eigen::VectorXd::Zero(0);
        const size_t targetFunctionRowIndex = isTargetFunctionModified()? simplexTable.rows() - 2 : simplexTable.rows() - 1;
        const size_t targetFunctionColIndex = simplexTable.cols() - 1;
        const double targetFunctionValue = simplexTable(targetFunctionRowIndex, targetFunctionColIndex);


        SimplexOutput output(
            SimplexResult(simplexSolution, iteration, targetFunctionValue, problemType), 
            input,
            isFeasible
        );
        output.writeResult();

        return output;
    }

    inline bool isTargetFunctionModified() const noexcept {
        return virtualVariableIndices.size() != 0;
    }
    
public:
    Simplex(const SimplexInput& simplexInput) : input(std::make_shared<SimplexInput>(simplexInput)) {}

    Simplex(SimplexInput&& simplexInput) : input(std::make_shared<SimplexInput>(std::move(simplexInput))) {}

    ~Simplex() = default;

    SimplexOutput solve(ProblemType problemType) {
        this->problemType = problemType;

        buildSimplexTable();

        excludeVirtualVariables();

        int iteration = 0;
        while (++iteration != SIMPLEX_MAX_ITERATIONS) {
            #ifdef __DEBUG__
                std::cout << "Iteration: " << iteration << '\n';
                std::cout << "Simplex table:\n" << simplexTable << '\n';
            #endif
            
            SimplexStepCode code = simplexStep();

            #ifdef __DEBUG__
                std::cout << "Simplex step state: ";
                switch(code) {
                    case SimplexStepCode::NEXT_STEP:
                        std::cout << "success.";
                        break;
                    case SimplexStepCode::NO_MAIN_COLUMN:
                        std::cout << "main column not found.";
                        break;
                    case SimplexStepCode::NO_MAIN_ROW:
                        std::cout << "main row not found.";
                        break;
                }
                std::cout << '\n';
            #endif

            if (code != SimplexStepCode::NEXT_STEP) {
                break;
            }
        }

        return generateOutput(iteration);
    }
};