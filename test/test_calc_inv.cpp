#include "doctest.h"

#include <Eigen/Dense>

// Local reimplementation of the update math from calc_R_inv_forward (non-GCTA).
static Eigen::MatrixXd update_inv_forward(
    const Eigen::MatrixXd& R_inv_pre,
    const Eigen::RowVectorXd& r_temp_row)
{
    Eigen::RowVectorXd temp_vector = r_temp_row * R_inv_pre;
    double temp_element = 1.0 / (1.0 - r_temp_row.dot(temp_vector));

    const int dim = static_cast<int>(R_inv_pre.rows());
    Eigen::MatrixXd R_inv_post = R_inv_pre + temp_element * temp_vector.transpose() * temp_vector;
    R_inv_post.conservativeResize(dim + 1, dim + 1);
    R_inv_post.topRightCorner(dim, 1) = -temp_element * temp_vector.transpose();
    R_inv_post.bottomLeftCorner(1, dim) = -temp_element * temp_vector;
    R_inv_post(dim, dim) = temp_element;

    return R_inv_post;
}

// Local reimplementation of the update math from calc_R_inv_backward (non-GCTA).
static Eigen::MatrixXd update_inv_backward(
    const Eigen::MatrixXd& R_inv_pre,
    int remove_index)
{
    Eigen::MatrixXd R_inv_post =
        R_inv_pre - R_inv_pre.col(remove_index) * R_inv_pre.row(remove_index) / R_inv_pre(remove_index, remove_index);
    const int dim = static_cast<int>(R_inv_post.rows());
    R_inv_post.block(remove_index, 0, dim - remove_index - 1, dim) =
        R_inv_post.block(remove_index + 1, 0, dim - remove_index - 1, dim);
    R_inv_post.conservativeResize(dim - 1, Eigen::NoChange);
    R_inv_post.block(0, remove_index, dim - 1, dim - remove_index - 1) =
        R_inv_post.block(0, remove_index + 1, dim - 1, dim - remove_index - 1);
    R_inv_post.conservativeResize(Eigen::NoChange, dim - 1);
    return R_inv_post;
}

TEST_CASE("forward inverse update matches direct inverse") {
    Eigen::MatrixXd R_pre(2, 2);
    R_pre << 1.0, 0.2,
             0.2, 1.0;

    Eigen::RowVectorXd r(2);
    r << 0.3, -0.1;

    Eigen::MatrixXd R_post(3, 3);
    R_post.topLeftCorner(2, 2) = R_pre;
    R_post.topRightCorner(2, 1) = r.transpose();
    R_post.bottomLeftCorner(1, 2) = r;
    R_post(2, 2) = 1.0;

    Eigen::MatrixXd R_inv_pre = R_pre.ldlt().solve(Eigen::MatrixXd::Identity(2, 2));
    Eigen::MatrixXd R_inv_expected = R_post.ldlt().solve(Eigen::MatrixXd::Identity(3, 3));
    Eigen::MatrixXd R_inv_post = update_inv_forward(R_inv_pre, r);

    CHECK(R_inv_post.isApprox(R_inv_expected, 1e-12));
}

TEST_CASE("backward inverse update matches direct inverse") {
    Eigen::MatrixXd R(3, 3);
    R << 1.0, 0.2, -0.1,
         0.2, 1.0, 0.3,
        -0.1, 0.3, 1.0;

    Eigen::MatrixXd R_inv = R.ldlt().solve(Eigen::MatrixXd::Identity(3, 3));

    const int remove_index = 1;
    Eigen::MatrixXd R_reduced(2, 2);
    R_reduced << 1.0, -0.1,
                 -0.1, 1.0;

    Eigen::MatrixXd R_inv_expected =
        R_reduced.ldlt().solve(Eigen::MatrixXd::Identity(2, 2));
    Eigen::MatrixXd R_inv_post = update_inv_backward(R_inv, remove_index);

    CHECK(R_inv_post.isApprox(R_inv_expected, 1e-12));
}
