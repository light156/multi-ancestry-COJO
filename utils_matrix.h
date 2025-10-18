#pragma once
#include <Eigen/Dense>
#include <vector>
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::PlainObjectBase;
using Eigen::DenseBase;
using Eigen::LDLT;


double median(const std::vector<double> &v);
double median(const ArrayXd &eigen_vector);
double calc_inner_product(const ArrayXd &vec1, const ArrayXd &vec2, bool if_keep_NA);

/*
// Fast Schur-complement-based block inverse update
bool calc_R_inverse_fast(
    const MatrixXd& R_inv_pre,
    const VectorXd& r_temp_vec,
    double lower_right_corner,
    double iter_colinear_threshold,
    MatrixXd& R_inv_post);
*/

// LDLT-based block inverse (GCTA-style)
bool calc_R_inverse_forward(
    const MatrixXd& R_pre,
    const VectorXd& r_temp_vec,
    double lower_right_corner,
    double iter_colinear_threshold,
    MatrixXd& R_post,
    MatrixXd& R_inv_post,
    bool do_check,
    const ArrayXd& scaling_vector=ArrayXd());

void calc_R_inverse_backward(
    const MatrixXd& R_pre,
    const MatrixXd& R_inv_pre,
    int remove_index,
    MatrixXd& R_post,
    MatrixXd& R_inv_post);

// template functions for matrix manipulation   
// Append a row: row must be 1 × N and match matrix.cols()
template <typename MatDerived, typename RowDerived>
inline void append_row(PlainObjectBase<MatDerived>& matrix, const DenseBase<RowDerived>& row_vec)
{
    if (row_vec.rows() != 1)
        throw std::invalid_argument("append_row: row must have 1 row (row vector).");
    if (matrix.rows() != 0 && matrix.cols() != row_vec.cols())
        throw std::invalid_argument("append_row: column count mismatch.");

    typename MatDerived::Index numRows = matrix.rows();

    matrix.conservativeResize(numRows + 1, row_vec.cols());
    matrix.row(numRows) = row_vec.derived();
}

// Append a column: col must be M × 1 and match matrix.rows()
template <typename MatDerived, typename ColDerived>
inline void append_column(PlainObjectBase<MatDerived>& matrix, const DenseBase<ColDerived>& col_vec)
{
    if (col_vec.cols() != 1)
        throw std::invalid_argument("append_column: col must have 1 column (column vector).");
    if (matrix.cols() != 0 && matrix.rows() != col_vec.rows())
        throw std::invalid_argument("append_column: row count mismatch.");

    typename MatDerived::Index numCols = matrix.cols();

    matrix.conservativeResize(col_vec.rows(), numCols + 1);
    matrix.col(numCols) = col_vec.derived();
}

// Remove a row by index (-1 for last row)
template <typename Derived>
inline void remove_row(PlainObjectBase<Derived> &matrix, typename Derived::Index index=-1) 
{
    typename Derived::Index numRows = matrix.rows() - 1;

    if (numRows <= 0 || matrix.cols() == 0) {
        matrix.resize(0, 0);
        return;
    }

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows - index) = matrix.bottomRows(numRows - index).eval();

    matrix.conservativeResize(numRows, Eigen::NoChange);
}

template <typename Derived>
inline void remove_column(PlainObjectBase<Derived> &matrix, typename Derived::Index index=-1) 
{
    typename Derived::Index numCols = matrix.cols() - 1;

    if (numCols <= 0 || matrix.rows() == 0) {
        matrix.resize(0, 0);
        return;
    }

    if (index != -1 && index < numCols)
        matrix.middleCols(index, numCols - index) = matrix.rightCols(numCols - index).eval();

    matrix.conservativeResize(Eigen::NoChange, numCols);
}
