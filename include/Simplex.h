#include <Eigen/Dense>
#include <iostream>

#include "SimplexResult.h"
#include "common.h"

enum class CompareSign {LESS_EQUAL, EQUAL, GREATER_EQUAL};
enum class SimplexStepCode {NEXT_STEP, NO_MAIN_COLUMN, NO_MAIN_ROW};

class Simplex final {
private:
    static constexpr int    SIMPLEX_MAX_ITERATIONS = 50;
    static constexpr double SIMPLEX_ACCURACY       = 1e-6;

    Eigen::VectorXi         compareSigns;
    Eigen::VectorXi         basicVariableIndices;
    Eigen::VectorXi         naturalVariableIndices;
    Eigen::VectorXi         virtualVariableIndices;

    Eigen::VectorXd         pricesVector;
    Eigen::VectorXd         boundsVector;

    Eigen::MatrixXd         simplexTable;
    Eigen::MatrixXd         boundsMatrix;

    ProblemType             problemType;

    void buildSimplexTable() {
        makeCoefficientsPositive();

        allocateMemory();

        findAllVariables();

        addConstraintsCoefficients();

        addTargetFunctionRow();

        if (!isTargetFunctionModified()) return;

        addArtificialTargetFunctionRow();
    }

    void makeCoefficientsPositive() {
        for (size_t row = 0; row < boundsMatrix.rows(); ++row) {
            if (boundsVector[row] >= 0) continue;
            boundsVector[row] *= -1.0;
            simplexTable.row(row) *= -1.0;
            
            if (CompareSign(compareSigns[row]) == CompareSign::EQUAL) continue;
            compareSigns[row] = static_cast<int>(CompareSign(compareSigns[row]) == CompareSign::LESS_EQUAL? CompareSign::GREATER_EQUAL : CompareSign::LESS_EQUAL);
        }

        #ifdef __DEBUG__
            std::cout << "Made coefficients positive\n";
        #endif
    }

    void allocateMemory() {
        size_t naturalNumber = pricesVector.rows() + boundsVector.rows(), basicNumber = 0, virtualNumber = 0;
        CompareSign compareSign;

        for (int value : compareSigns) {
            compareSign = static_cast<CompareSign>(value);
            if (compareSign == CompareSign::LESS_EQUAL) {
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

        simplexTable = Eigen::MatrixXd::Zero(boundsVector.rows() + 1 + isTargetFunctionModified(), naturalNumber + virtualNumber + 1);

        #ifdef __DEBUG__
            std::cout << "Allocated simplex table\n";
            std::cout << "Rows: " << boundsVector.rows() + 1 + isTargetFunctionModified() << "\n";
            std::cout << "Cols: " << naturalNumber + virtualNumber + 1 << "\n";
        #endif
    }

    void findAllVariables() {
        size_t naturalVariablesCount = 0;
        size_t basicVariablesCount = 0;
        size_t virtualVariablesCount = 0;

        for (size_t index = 0; index < pricesVector.rows(); ++index) {
            naturalVariableIndices[naturalVariablesCount++] = index;
        }

        #ifdef __DEBUG__
            std::cout << "Added variables from target function\n";
        #endif

        size_t lastFilledColIndex = boundsMatrix.cols() - 1;

        int additionalVariableIndex, virtualVariableIndex;

        for (size_t index = 0; index < boundsMatrix.rows(); ++index) {
            createArtificialVariables(index, additionalVariableIndex, virtualVariableIndex, lastFilledColIndex);

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

    void createArtificialVariables(int inequalityIndex, 
                                int& additionalVariableIndex, 
                                int& virtualVariableIndex,
                                size_t& lastFilledColIndex) {
        if (CompareSign(compareSigns[inequalityIndex]) == CompareSign::GREATER_EQUAL) {
            simplexTable(inequalityIndex, ++lastFilledColIndex) = -1.0;
            simplexTable(inequalityIndex, ++lastFilledColIndex) = 1.0;
            additionalVariableIndex = lastFilledColIndex - 1;
            virtualVariableIndex = lastFilledColIndex;
        } else {
            simplexTable(inequalityIndex, ++lastFilledColIndex) = 1.0;
            additionalVariableIndex = lastFilledColIndex;
            virtualVariableIndex = CompareSign(compareSigns[inequalityIndex]) == CompareSign::EQUAL? additionalVariableIndex : -1;
        }
    }

    void addConstraintsCoefficients() {
        for (size_t row = 0; row < boundsMatrix.rows(); ++row) {
            for (size_t col = 0; col < boundsMatrix.cols(); ++col) {
                simplexTable(row, col) = boundsMatrix(row, col);
            }
        }

        #ifdef __DEBUG__
            std::cout << "Carried values from bounds matrix to simplex table\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif

        simplexTable.block(0, simplexTable.cols() - 1, boundsVector.rows(), 1) = boundsVector;

        #ifdef __DEBUG__
            std::cout << "Carried free coefficients from bounds vector to simplex table\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif
    }

    void addTargetFunctionRow() {
        Eigen::VectorXd targetFunctionRow = Eigen::VectorXd::Zero(simplexTable.cols());
        
        for (size_t index = 0; index < pricesVector.rows(); ++index) {
            targetFunctionRow[index] = problemType == ProblemType::MAX ? -pricesVector[index] : pricesVector[index];
        }

        simplexTable.row(boundsMatrix.rows()) = targetFunctionRow;

        #ifdef __DEBUG__
            std::cout << "Adding target function row\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif
    }

    void addArtificialTargetFunctionRow() {
        Eigen::VectorXd newTargetFunctionRow = Eigen::VectorXd::Zero(simplexTable.cols());
        
        for (size_t virtualVariableIndex : virtualVariableIndices) {
            newTargetFunctionRow[virtualVariableIndex] = 1.0;
        }

        simplexTable.row(boundsMatrix.rows() + 1) = newTargetFunctionRow;

        #ifdef __DEBUG__
            std::cout << "Added artificial target function row\n";
            std::cout << "Simplex table:\n" << simplexTable << '\n';
        #endif
    }

    void excludeVirtualVariables() {
        if (!isTargetFunctionModified()) return;

        const size_t lastRowId = simplexTable.rows() - 1;

        for (size_t virtualVariableIndex : virtualVariableIndices) {
            for (size_t row = 0; row < simplexTable.rows(); ++row) {
                if (simplexTable(row, virtualVariableIndex) == 0) continue;
                double arg = simplexTable(lastRowId, virtualVariableIndex) / simplexTable(row, virtualVariableIndex);
                simplexTable.row(lastRowId) -= arg * simplexTable.row(row);
                break;
            }
        }

        #ifdef __DEBUG__
            std::cout << "Excluded virtual variables\n";
        #endif
    }

    bool isPlanOptimal() {
        bool isOptimal = true;
        for (int value : simplexTable.row(simplexTable.rows() - 1)) {
            isOptimal = value < 0;
            if (!isOptimal) break;
        }

        if (!isTargetFunctionModified()) return isOptimal;
        
        for (int value : simplexTable.row(simplexTable.rows() - 2)) {
            isOptimal &= value < 0;
            if (!isOptimal) break;
        }

        return isOptimal;
    }

    int getMainColumn() {
        Eigen::RowVectorXd row = simplexTable.row(simplexTable.rows() - 1);
        double minimumNotPositive = 0;
        int mainCol = -1;

        for (size_t index = 0; index < row.cols() - 1; ++index) {
            if (row[index] >= minimumNotPositive) continue;
            minimumNotPositive = row[index];
            mainCol = index;
        }

        if (!isTargetFunctionModified() || mainCol != -1) return mainCol;

        row = simplexTable.row(simplexTable.rows() - 2);

        for (size_t index : naturalVariableIndices) {
            if (row[index] >= minimumNotPositive) continue;
            minimumNotPositive = row[index];
            mainCol = index;
        }

        return mainCol;
    }

    int getMainRow(int mainCol) {
        double minimumValue = std::numeric_limits<double>::max(), mainElement;
        int mainRow = -1, lastColIndex = simplexTable.cols() - 1, targetFunctionRowIndex = isTargetFunctionModified()? simplexTable.rows() - 2 : simplexTable.rows() - 1;

        for (size_t index = 0; index < targetFunctionRowIndex; ++index) {
            mainElement = simplexTable(index, mainCol);
            if (mainElement < 0) continue;
            if (simplexTable(index, lastColIndex) / mainElement > minimumValue) continue;

            minimumValue = simplexTable(index, lastColIndex) / mainElement;
            mainRow = index; 
        }

        return mainRow;
    }

    bool validateSolution() {
        double targetFunctionValue = 0;

        const size_t targetFunctionRowIndex = isTargetFunctionModified()? simplexTable.rows() - 2 : simplexTable.rows() - 1;
        const size_t targetFunctionColIndex = simplexTable.cols() - 1;

        for (size_t index = 0; index < basicVariableIndices.rows(); ++index) {
            if (basicVariableIndices[index] >= pricesVector.rows()) {
                continue;
            }

            targetFunctionValue += simplexTable(index, targetFunctionColIndex) * pricesVector[basicVariableIndices[index]];
        }

        double sign = problemType == ProblemType::MAX? -1.0 : 1.0;
        double targetFunctionTableValue = simplexTable.row(targetFunctionRowIndex)[targetFunctionColIndex];
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

    Eigen::VectorXd getSimplexSolution(bool onlyNaturalArgs) {
        Eigen::VectorXd solution = Eigen::VectorXd::Zero(onlyNaturalArgs? pricesVector.rows() : simplexTable.cols() - 1);

        for (size_t index = 0; index < basicVariableIndices.rows(); ++index) {
            if (basicVariableIndices[index] >= solution.rows()) continue;
            solution[basicVariableIndices[index]] = simplexTable(index, simplexTable.cols() - 1);
        }

        return solution;
    }

    SimplexStepCode simplexStep() {
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

    SimplexResult generateResult(const int iteration) {
        Eigen::VectorXd simplexSolution = validateSolution()? getSimplexSolution(true) : Eigen::VectorXd::Zero(0);
        const size_t targetFunctionRowIndex = isTargetFunctionModified()? simplexTable.rows() - 2 : simplexTable.rows() - 1;
        const size_t targetFunctionColIndex = simplexTable.cols() - 1;
        const double targetFunctionValue = simplexTable.row(targetFunctionRowIndex)[targetFunctionColIndex];

        SimplexResult result(simplexSolution, targetFunctionValue, this->problemType, iteration);
        result.writeResult(compareSigns, pricesVector, boundsVector, boundsMatrix);

        return result;
    }

    inline bool isTargetFunctionModified() {
        return virtualVariableIndices.size() != 0;
    }
    
public:
    Simplex(const Eigen::MatrixXd& boundsMatrix, 
            const Eigen::VectorXd& boundsVector, 
            const Eigen::VectorXd& pricesVector, 
            const Eigen::VectorXi& compareSigns) {
        if (boundsMatrix.rows() != boundsVector.rows()) {
            throw new std::invalid_argument("Error, when creating an instance of Simplex: number of rows in boundsMatrix is not equal to number of elements in boundsVector");
        }

        if (boundsMatrix.cols() != pricesVector.rows()) {
            throw new std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in pricesVector");
        }

        if (boundsVector.rows() != compareSigns.rows()) {
            throw new std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in compareSigns");
        }

        this->boundsMatrix = boundsMatrix;
        this->boundsVector = boundsVector;
        this->pricesVector = pricesVector;
        this->compareSigns = compareSigns;
    }

    Simplex(Eigen::MatrixXd&& boundsMatrix, 
            Eigen::VectorXd&& boundsVector, 
            Eigen::VectorXd&& pricesVector, 
            Eigen::VectorXi&& compareSigns) {
        if (boundsMatrix.rows() != boundsVector.rows()) {
            throw new std::invalid_argument("Error, when creating an instance of Simplex: number of rows in boundsMatrix is not equal to number of elements in boundsVector");
        }

        if (boundsMatrix.cols() != pricesVector.rows()) {
            throw new std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in pricesVector");
        }

        if (boundsVector.rows() != compareSigns.rows()) {
            throw new std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in compareSigns");
        }

        this->boundsMatrix = std::move(boundsMatrix);
        this->boundsVector = std::move(boundsVector);
        this->pricesVector = std::move(pricesVector);
        this->compareSigns = std::move(compareSigns);
    }

    ~Simplex() = default;

    SimplexResult solve(ProblemType problemType) {
        this->problemType = problemType;

        buildSimplexTable();

        excludeVirtualVariables();

        int iteration = 0;
        while ((!isPlanOptimal()) && (iteration != SIMPLEX_MAX_ITERATIONS)) {
            ++iteration;
            #ifdef __DEBUG__
                std::cout << "Iteration: " << iteration << '\n';
                std::cout << "Simplex table:\n" << simplexTable << '\n';
            #endif
            
            if(simplexStep() != SimplexStepCode::NEXT_STEP) {
                break;
            }
        }

        return generateResult(iteration);
    }

    const Eigen::VectorXi& getCompareSigns() {
        return compareSigns;
    }

    const Eigen::VectorXd& getPricesVector() {
        return pricesVector;
    }

    const Eigen::VectorXd& getBoundsVector() {
        return boundsVector;
    }

    const Eigen::MatrixXd& getBoundsMatrix() {
        return boundsMatrix;
    }
};