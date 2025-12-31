#pragma once

#include <external/Eigen/Dense>
#include "external/fast_float.h"
#include <vector>
#include <cmath>
#include <algorithm>

using Eigen::PlainObjectBase;
using Eigen::DenseBase;
using std::invalid_argument;
using std::vector;
using std::string;
using std::pair;


// template functions for matrix manipulation   
// Append a row: row must be 1 × N and match matrix.cols()
template <typename MatDerived, typename RowDerived>
inline void append_row(PlainObjectBase<MatDerived>& matrix, const DenseBase<RowDerived>& row_vec)
{
    if (row_vec.rows() != 1)
        throw invalid_argument("append_row: row must have 1 row (row vector).");
    if (matrix.rows() != 0 && matrix.cols() != row_vec.cols())
        throw invalid_argument("append_row: column count mismatch.");

    typename MatDerived::Index numRows = matrix.rows();

    matrix.conservativeResize(numRows + 1, row_vec.cols());
    matrix.row(numRows) = row_vec.derived();
}


// Append a column: col must be M × 1 and match matrix.rows()
template <typename MatDerived, typename ColDerived>
inline void append_column(PlainObjectBase<MatDerived>& matrix, const DenseBase<ColDerived>& col_vec)
{
    if (col_vec.cols() != 1)
        throw invalid_argument("append_column: col must have 1 column (column vector).");
    if (matrix.cols() != 0 && matrix.rows() != col_vec.rows())
        throw invalid_argument("append_column: row count mismatch.");

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


inline double median(std::vector<double>& v) {
    const size_t n = v.size();
    if (n == 0) 
        throw invalid_argument("median: empty vector");

    const size_t mid = n / 2;

    // First nth_element: put the upper median in place
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double upper = v[mid];

    if (n % 2 == 1) {
        return upper; // odd size, single middle value
    } else {
        std::nth_element(v.begin(), v.begin() + mid - 1, v.end()); // even size, need lower median too
        double lower = v[mid - 1];
        return 0.5 * (lower + upper); 
    }
}


inline double median(const Eigen::ArrayXd &eigen_vector)
{
    vector<double> v(eigen_vector.data(), eigen_vector.data() + eigen_vector.size());
    return median(v);
}


inline void skip_delim(const char*& pt) {
    while (*pt && (*pt == ' ' || *pt == '\t' || *pt == '\r')) pt++;
}


inline void skip_token(const char*& pt) {
    while (*pt && *pt != ' ' && *pt != '\t' && *pt != '\r') pt++;
}


inline void parse_string(const char*& pt, string& out, bool to_upper=false) 
{
    skip_delim(pt);
    const char* start = pt;
    
    skip_token(pt);
    size_t len = pt - start;
    
    out.resize(len);
    std::memcpy(&out[0], start, len);

    if (to_upper) {
        for (size_t i = 0; i < len; i++) {
            char c = out[i];
            out[i] = (c >= 'a' && c <= 'z') ? (c - 32) : c;
        }
    }
}


template<typename T>
inline bool parse_num(const char*& pt, T& out) {
    skip_delim(pt);
    const char* start = pt;

    skip_token(pt);
    auto answer = fast_float::from_chars(start, pt, out);

    if (answer.ec != std::errc() || answer.ptr != pt) return false;
    if (std::isnan(out) || std::isinf(out)) return false;
    return true;
}


inline vector<pair<string,int>>::iterator fast_lookup(vector<pair<string,int>>& table, const string& key)
{
    auto iter = std::lower_bound(table.begin(), table.end(), key,
        [](const pair<string,int>& a, const string& key) { return a.first < key;});

    if (iter == table.end() || iter->first != key || iter->second == -1) return table.end();
    return iter;
}
