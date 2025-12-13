#pragma once

#include <Eigen/Dense>
#include "common.h"

class SimplexInput final {
private:
    Eigen::VectorXi compareSigns;
    Eigen::VectorXd pricesVector;
    Eigen::VectorXd boundsVector;
    Eigen::MatrixXd boundsMatrix;
public:
    SimplexInput(const Eigen::MatrixXd& boundsMatrix, 
            const Eigen::VectorXd& boundsVector, 
            const Eigen::VectorXd& pricesVector, 
            const Eigen::VectorXi& compareSigns) {
        if (boundsMatrix.rows() != boundsVector.rows()) {
            throw std::invalid_argument("Error, when creating an instance of Simplex: number of rows in boundsMatrix is not equal to number of elements in boundsVector");
        }

        if (boundsMatrix.cols() != pricesVector.rows()) {
            throw std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in pricesVector");
        }

        if (boundsVector.rows() != compareSigns.rows()) {
            throw std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in compareSigns");
        }

        this->boundsMatrix = boundsMatrix;
        this->boundsVector = boundsVector;
        this->pricesVector = pricesVector;
        this->compareSigns = compareSigns;
    }

    SimplexInput(Eigen::MatrixXd&& boundsMatrix, 
            Eigen::VectorXd&& boundsVector, 
            Eigen::VectorXd&& pricesVector, 
            Eigen::VectorXi&& compareSigns) {
        if (boundsMatrix.rows() != boundsVector.rows()) {
            throw std::invalid_argument("Error, when creating an instance of Simplex: number of rows in boundsMatrix is not equal to number of elements in boundsVector");
        }

        if (boundsMatrix.cols() != pricesVector.rows()) {
            throw std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in pricesVector");
        }

        if (boundsVector.rows() != compareSigns.rows()) {
            throw std::invalid_argument("Error, when creating an instance of Simplex: number of cols in boundsMatrix is not equal to number of elements in compareSigns");
        }

        this->boundsMatrix = std::move(boundsMatrix);
        this->boundsVector = std::move(boundsVector);
        this->pricesVector = std::move(pricesVector);
        this->compareSigns = std::move(compareSigns);
    }

    SimplexInput(const SimplexInput& simplexInput) noexcept
    : compareSigns(simplexInput.compareSigns),
    pricesVector(simplexInput.pricesVector),
    boundsVector(simplexInput.boundsVector),
    boundsMatrix(simplexInput.boundsMatrix) {}

    SimplexInput(SimplexInput&& simplexInput) noexcept 
    : compareSigns(std::move(simplexInput.compareSigns)),
    pricesVector(std::move(simplexInput.pricesVector)),
    boundsVector(std::move(simplexInput.boundsVector)),
    boundsMatrix(std::move(simplexInput.boundsMatrix)) {}

    ~SimplexInput() = default;

    SimplexInput& operator=(const SimplexInput& other) {
        if (this != &other) {
            compareSigns = other.compareSigns;
            pricesVector = other.pricesVector;
            boundsVector = other.boundsVector;
            boundsMatrix = other.boundsMatrix;
        }

        return *this;
    }

    inline const Eigen::VectorXi& getCompareSigns() const noexcept {
        return compareSigns;
    }

    inline const Eigen::VectorXd& getPricesVector() const noexcept {
        return pricesVector;
    }

    inline const Eigen::VectorXd& getBoundsVector() const noexcept {
        return boundsVector;
    }

    inline const Eigen::MatrixXd& getBoundsMatrix() const noexcept {
        return boundsMatrix;
    }

    friend class Simplex;
    friend class SimplexOutput;
};