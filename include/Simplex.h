#include <Eigen/Dense>
#include <iostream>

enum class CompareSign {LESS_EQUAL, EQUAL, GREATER_EQUAL};
enum class ProblemType {MAX, MIN};

class Simplex {
private:
    static constexpr int SIMPLEX_MAX_ITERATIONS = 50;
    static constexpr double SIMPLEX_ACCURACY       = 1e-6;

    Eigen::VectorXi compareSigns;
    Eigen::VectorXd pricesVector;
    Eigen::VectorXd boundsVector;
    Eigen::VectorXi basicVariableIndices;
    Eigen::VectorXi naturalVariableIndices;
    Eigen::VectorXi virtualVariableIndices;

    Eigen::MatrixXd simplexTable;
    Eigen::MatrixXd boundsMatrix;

    ProblemType     problemType = ProblemType::MAX;

    bool isTargetFunctionModified() {
        return virtualVariableIndices.size() != 0;
    }

    void buildSimplexTable() {
        #ifdef __DEBUG__
            std::cout << "Making coefficients positive\n";
        #endif
        for (size_t row = 0; row < simplexTable.rows(); ++row) {
            if (boundsVector[row] >= 0) continue;
            boundsVector[row] *= -1.0;
            simplexTable.row(row) *= -1.0;

            if (CompareSign(compareSigns[row]) == CompareSign::EQUAL) continue;
            compareSigns[row] = static_cast<int>(CompareSign(compareSigns[row]) == CompareSign::LESS_EQUAL? CompareSign::GREATER_EQUAL : CompareSign::LESS_EQUAL);
        }

        #ifdef __DEBUG__
            std::cout << "Precompute all dimensions of vectors and matrix\n";
        #endif
        allocateMemory();

        #ifdef __DEBUG__
            std::cout << "Taking values from bounds matrix to simplex table\n";
        #endif
        for (size_t row = 0; row < boundsMatrix.rows(); ++row) {
            for (size_t col = 0; col < boundsMatrix.cols(); ++col) {
                simplexTable(row, col) = boundsMatrix(row, col);
            }
        }

        size_t naturalVariablesCount = 0;
        size_t basicVariablesCount = 0;
        size_t virtualVariablesCount = 0;
        #ifdef __DEBUG__
            std::cout << "Adding natural variables\n";
        #endif
        for (size_t index = 0; index < pricesVector.rows(); ++index) {
            naturalVariableIndices[naturalVariablesCount++] = index;
        }

        size_t lastFilledRowIndex = boundsMatrix.rows() - 1;
        size_t lastFilledColIndex = boundsMatrix.cols() - 1;

        int additionalVariableIndex, virtualVariableIndex;

        #ifdef __DEBUG__
            std::cout << "Creating new variables\n";
        #endif
        for (size_t index = 0; index < compareSigns.rows(); ++index) {
            createArtificialVariables(index, additionalVariableIndex, virtualVariableIndex, lastFilledColIndex);

            naturalVariableIndices[naturalVariablesCount++] = additionalVariableIndex;

            if (virtualVariableIndex != -1) {
                basicVariableIndices[basicVariablesCount++] = virtualVariableIndex;
                virtualVariableIndices[virtualVariablesCount++] = virtualVariableIndex;
            }

            basicVariableIndices[basicVariablesCount++] = additionalVariableIndex;
        }

        #ifdef __DEBUG__
            std::cout << "Adding free coefficients\n";
        #endif
        simplexTable.block(0, ++lastFilledColIndex, boundsVector.rows(), 1) = boundsVector;
        
        #ifdef __DEBUG__
            std::cout << "Adding target function row\n";
        #endif
        Eigen::VectorXd targetFunctionRow = Eigen::VectorXd::Zero(simplexTable.cols());
        for (size_t index = 0; index < pricesVector.rows(); ++index) {
            targetFunctionRow[index] = problemType == ProblemType::MAX? -pricesVector[index] : pricesVector[index];
        }

        simplexTable.row(++lastFilledRowIndex) = targetFunctionRow;

        if (!isTargetFunctionModified()) return;

        #ifdef __DEBUG__
            std::cout << "Adding new target function row for slacks\n";
        #endif
        Eigen::VectorXd newTargetFunctionRow = Eigen::VectorXd::Zero(simplexTable.cols());
        for (size_t virtualVariableIndex : virtualVariableIndices) {
            newTargetFunctionRow[virtualVariableIndex] = 1.0;
        }

        simplexTable.row(++lastFilledRowIndex) = newTargetFunctionRow;
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

        #ifdef __DEBUG__
            std::cout << "Allocating vectors\n";
            std::cout << "Natural variables amount: " << naturalNumber << "\n";
            std::cout << "Basic variables amount: " << basicNumber << "\n";
            std::cout << "Virtual variables amount: " << virtualNumber << "\n";
        #endif
        basicVariableIndices = Eigen::VectorXi(basicNumber);
        naturalVariableIndices = Eigen::VectorXi(naturalNumber);
        virtualVariableIndices = Eigen::VectorXi(virtualNumber);
        #ifdef __DEBUG__
            std::cout << "Allocating simplex table\n";
            std::cout << "Rows: " << boundsVector.rows() + 1 + isTargetFunctionModified() << "\n";
            std::cout << "Cols: " << naturalNumber + virtualNumber + 1 << "\n";
        #endif
        simplexTable = Eigen::MatrixXd::Zero(boundsVector.rows() + 1 + isTargetFunctionModified(), naturalNumber + virtualNumber + 1);
    }

    void createArtificialVariables(int inequalityIndex, 
                                int& additionalVariableIndex, 
                                int& virtualVariableIndex,
                                size_t& lastFilledColIndex) {
        if (CompareSign(compareSigns[inequalityIndex]) == CompareSign::GREATER_EQUAL) {
            simplexTable.col(++lastFilledColIndex) = Eigen::VectorXd::Zero(simplexTable.rows());
            simplexTable.col(++lastFilledColIndex) = Eigen::VectorXd::Zero(simplexTable.rows());
            simplexTable(inequalityIndex, lastFilledColIndex) = 1.0;
            simplexTable(inequalityIndex, lastFilledColIndex - 1) = -1.0;
            additionalVariableIndex = lastFilledColIndex - 1;
            virtualVariableIndex = lastFilledColIndex;
        } else {
            simplexTable.col(++lastFilledColIndex) = Eigen::VectorXd::Zero(simplexTable.rows());
            simplexTable(inequalityIndex, lastFilledColIndex) = 1.0;
            additionalVariableIndex = lastFilledColIndex;
            virtualVariableIndex = CompareSign(compareSigns[inequalityIndex]) == CompareSign::EQUAL? additionalVariableIndex : -1;
        }
    }

    bool excludeVirtualVariables() {
        if (!isTargetFunctionModified()) return false;

        const int lastRowId = simplexTable.rows() - 1;

        for (size_t virtualVariableIndex : virtualVariableIndices) {
            for (size_t row = 0; row < simplexTable.rows(); ++row) {
                if (simplexTable(row, virtualVariableIndex) == 0) continue;
                double arg = simplexTable(lastRowId, virtualVariableIndex) / simplexTable(row, virtualVariableIndex);
                simplexTable.row(lastRowId) -= arg * simplexTable.row(row);
                break;
            }
        }

        return true;
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

    int getMainCol() {
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

public:

    Simplex(const Eigen::MatrixXd& boundsMatrix, const Eigen::VectorXd& boundsVector, const Eigen::VectorXd& pricesVector, const Eigen::VectorXi& compareSigns) {
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

    Simplex(Eigen::MatrixXd&& boundsMatrix, Eigen::VectorXd&& boundsVector, Eigen::VectorXd&& pricesVector, Eigen::VectorXi&& compareSigns) {
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

    Eigen::VectorXd solve(ProblemType problemType) {
        this->problemType = problemType;

        #ifdef __DEBUG__
            std::cout << "Start building simplex table\n";
        #endif
        buildSimplexTable();
        #ifdef __DEBUG__
            std::cout << "End building simplex table\n";
        #endif

        double mainElement;
        int mainRow, mainCol, iteration = 0;
        Eigen::VectorXd solution (pricesVector.size());

        if (excludeVirtualVariables()) {
            std::cout << "Virtual variables was removed from simplex table\n";
        }

        #ifdef __DEBUG__
            std::cout << "Starting main algorithm\n";
        #endif
        while ((!isPlanOptimal()) || (iteration != SIMPLEX_MAX_ITERATIONS)) {
            ++iteration;
            #ifdef __DEBUG__
                std::cout << "Iteration: " << iteration << '\n';
                std::cout << "Simplex table:\n" << simplexTable << '\n';
                std::cout << "Getting main column\n";
            #endif
            mainCol = getMainCol();

            if (mainCol == -1) {
                break;
            }

            #ifdef __DEBUG__
                std::cout << "Getting main row\n";
            #endif
            mainRow = getMainRow(mainCol);

            if (mainRow == -1) {
                return solution;
            }

            basicVariableIndices[mainRow] = mainCol;

            mainElement = simplexTable(mainRow, mainCol);

            simplexTable.row(mainRow) /= mainElement;

            for (size_t index = 0; index < simplexTable.rows(); ++index) {
                if (index == mainRow) continue;
                simplexTable.row(index) -= simplexTable(index, mainCol) * simplexTable.row(mainRow);
            }
        }

        #ifdef __DEBUG__
            std::cout << "Starting validating solution\n";
        #endif
        if (validateSolution()) {
            #ifdef __DEBUG__
                std::cout << "Solution is valid. Getting the answer\n";
            #endif
            solution = getSimplexSolution(true);
        }
        return solution;
    }
};