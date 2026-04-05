#pragma once

#include <external/Eigen/Dense>
#include "external/fast_float.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using Eigen::PlainObjectBase;
using Eigen::DenseBase;
using std::vector;
using std::string;
using std::pair;


// template functions for matrix manipulation
// Insert a row at position index (-1 = append at end)
template <typename MatDerived, typename RowDerived>
inline void insert_row(PlainObjectBase<MatDerived>& matrix, const DenseBase<RowDerived>& row_vec,
                       typename MatDerived::Index index = -1)
{
    if (row_vec.rows() != 1)
        throw std::invalid_argument("insert_row: row must have 1 row (row vector).");
    if (matrix.rows() != 0 && matrix.cols() != row_vec.cols())
        throw std::invalid_argument("insert_row: column count mismatch.");

    typename MatDerived::Index numRows = matrix.rows();
    if (index == -1) index = numRows;
    if (index < 0 || index > numRows)
        throw std::invalid_argument("insert_row: index out of range.");

    matrix.conservativeResize(numRows + 1, row_vec.cols());
    if (index < numRows)
        matrix.bottomRows(numRows - index) = matrix.middleRows(index, numRows - index).eval();
    matrix.row(index) = row_vec.derived();
}


// Insert a column at position index (-1 = append at end)
template <typename MatDerived, typename ColDerived>
inline void insert_column(PlainObjectBase<MatDerived>& matrix, const DenseBase<ColDerived>& col_vec,
                          typename MatDerived::Index index = -1)
{
    if (col_vec.cols() != 1)
        throw std::invalid_argument("insert_column: col must have 1 column (column vector).");
    if (matrix.cols() != 0 && matrix.rows() != col_vec.rows())
        throw std::invalid_argument("insert_column: row count mismatch.");

    typename MatDerived::Index numCols = matrix.cols();
    if (index == -1) index = numCols;
    if (index < 0 || index > numCols)
        throw std::invalid_argument("insert_column: index out of range.");

    matrix.conservativeResize(col_vec.rows(), numCols + 1);
    if (index < numCols)
        matrix.rightCols(numCols - index) = matrix.middleCols(index, numCols - index).eval();
    matrix.col(index) = col_vec.derived();
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
        throw std::invalid_argument("median: empty vector");

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
    vector<double> v(eigen_vector.size());
    for (auto i = 0; i < eigen_vector.size(); i++)
        v[i] = eigen_vector[i];
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
