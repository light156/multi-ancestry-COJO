#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
using Eigen::PlainObjectBase;
using Eigen::DenseBase;


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


inline double median(const std::vector<double> &v)
{
    int size = v.size();
	if (size == 0) 
        throw std::invalid_argument("median: empty vector");

    std::vector<double> b(v);
    std::sort(b.begin(), b.end());
	return (size%2==1) ? b[(size-1)/2] : (b[size/2]+b[size/2-1])/2;
}


inline double median(const Eigen::ArrayXd &eigen_vector)
{
    std::vector<double> v(eigen_vector.data(), eigen_vector.data() + eigen_vector.size());
    return median(v);
}


inline void skip_delim(const char*& p) {
    while (*p && (*p == ' ' || *p == '\t')) p++;
}


inline void skip_token(const char*& p) {
    while (*p && *p != ' ' && *p != '\t') p++;
}


inline void parse_string(const char*& p, std::string& out, bool to_upper=false) {
    skip_delim(p);
    const char* start = p;
    
    skip_token(p);
    if (start == p) return;
    
    size_t len = p - start;
    out.resize(len);
    std::memcpy(&out[0], start, len);

    if (to_upper) {
        for (size_t i = 0; i < len; i++) {
            char c = out[i];
            out[i] = (c >= 'a' && c <= 'z') ? (c - 32) : c;
        }
    }
}


inline void parse_int(const char*& p, int& out) {
    char* end;
    out = strtol(p, &end, 10);
    p = end;
    skip_delim(p);
}


inline bool parse_double(const char*& p, double& out) {
    char* end;
    out = strtod(p, &end);
    if (end == p) return false; // if end == p, it's NA, ".", end of line, junk, whatever, just consider it invalid
    p = end;
    return true;
}
