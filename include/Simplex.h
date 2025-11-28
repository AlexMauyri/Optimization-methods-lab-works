#include <Eigen/Dense>

enum class CompareSign {LESS_EQUAL=-1, EQUAL, GREATER_EQUAL};
enum class ProblemType {MAX, MIN};

class Simplex {
private:
    Eigen::VectorXi compare_signs;
    Eigen::VectorXd prices_vector;
    Eigen::VectorXd bounds_vector;
    Eigen::VectorXd indices_of_basic_variables;
    Eigen::VectorXd indices_of_free_variables;
    Eigen::MatrixXd simplex_table;
    Eigen::MatrixXd bounds_matrix;

    ProblemType     problem_type = ProblemType::MAX;
public:
    Simplex(const Eigen::MatrixXd& bounds_matrix, const Eigen::VectorXd& bounds_vector, const Eigen::VectorXd& prices_vector, const Eigen::VectorXi& compare_signs) {
        if (bounds_matrix.rows() != bounds_vector.cols()) {

        }

        if (bounds_matrix.cols() != prices_vector.cols()) {

        }

        if (bounds_vector.cols() != compare_signs.cols()) {

        }

        this->bounds_matrix = bounds_matrix;
        this->bounds_vector = bounds_vector;
        this->prices_vector = prices_vector;
        this->compare_signs = compare_signs;
    }

    Simplex(Eigen::MatrixXd&& bounds_matrix, Eigen::VectorXd&& bounds_vector, Eigen::VectorXd&& prices_vector, Eigen::VectorXi&& compare_signs);
};