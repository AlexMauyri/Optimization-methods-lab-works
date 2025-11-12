#pragma once

#include <functional>
#include <vector>
#include <Eigen/Dense>

#include "PenaltyFunction.h"
#include "multi_dim.h"
#include "search_result_nd.h"
#include "common.h"

enum class ConstraintType {EQUALITY, INEQUALITY};

class OuterPenaltyFunction final : public PenaltyFunction {
private:
    std::vector<function_nd> equalities;

    double computeTargetFunctionWithPenalty(const Eigen::VectorXd& start) const override { return target_function(start) + outerPenalty(start); }

    inline double outerPenalty(const Eigen::VectorXd& start) const { return equalityPenalty(start) + inequalityPenalty(start); }

    inline double equalityPenalty(const Eigen::VectorXd& start) const {
        return transformConstraints(start, ConstraintType::EQUALITY, 
            [this](double value) {
                return equalityValueTransform(value);
            }
        );
    }

    inline double inequalityPenalty(const Eigen::VectorXd& start) const {
        return transformConstraints(start, ConstraintType::INEQUALITY, 
            [this](double value) {
                return inequalityValueTransform(value);
            }
        );
    }

    double transformConstraints(const Eigen::VectorXd& start, 
                                ConstraintType type, 
                                const std::function<double(double)>& transform) const {
        const std::vector<function_nd>& constraints = (type == ConstraintType::EQUALITY)? this->equalities : this->inequalities;
        
        double accum = 0.0;
        for (const auto& constraint : constraints) {
            accum += transform(constraint(start));
        }
        return accum;
    }

    inline double equalityValueTransform(double value) const { return value * value; }

    inline double inequalityValueTransform(double value) const { 
        double not_negative = std::max(value, 0.0);
        return not_negative * not_negative;
    }

public:
    explicit OuterPenaltyFunction(const function_nd& target_function) noexcept 
        : PenaltyFunction(target_function) {}

    OuterPenaltyFunction(const OuterPenaltyFunction& other) noexcept 
        : PenaltyFunction(other)
        , equalities(other.getEqualities()) {}

    OuterPenaltyFunction(OuterPenaltyFunction&& other) noexcept 
        : PenaltyFunction(std::move(other))
        , equalities(other.getEqualities()) {}

    ~OuterPenaltyFunction() { equalities.clear(); }

    OuterPenaltyFunction& operator=(const OuterPenaltyFunction& other) noexcept {
        if (this != &other) {
            this->equalities = other.equalities;
            PenaltyFunction::operator=(other);
        }
        return *this;
    }

    OuterPenaltyFunction& operator=(OuterPenaltyFunction&& other) noexcept {
        if (this != &other) {
            this->equalities = other.equalities;
            PenaltyFunction::operator=(std::move(other));
        }
        return *this;
    }

    inline size_t getAmountOfEqualities() const { return equalities.size(); }

    inline const std::vector<function_nd>& getEqualities() const { return equalities; }

    inline void addEquality(function_nd equality) { equalities.push_back(equality); }

    inline void deleteEquality(size_t index) {
        if (index >= equalities.size()) {
            throw std::out_of_range("Given index is out of range.");
        }

        equalities.erase(equalities.begin() + index);
    }

    inline void clearEqualities() noexcept { equalities.clear(); }
};