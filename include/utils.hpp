#pragma once

#include <Eigen/Dense>
using Eigen::PlainObjectBase;
using Eigen::DenseBase;

#include <vector>
#include <sstream>
#include <cmath>
using std::vector;
using std::invalid_argument;


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


inline double median(const vector<double> &v)
{
    int size = v.size();
	if (size == 0) 
        throw invalid_argument("median: empty vector");

    vector<double> b(v);
    sort(b.begin(), b.end());
	return (size%2==1) ? b[(size-1)/2] : (b[size/2]+b[size/2-1])/2;
}


inline double median(const Eigen::ArrayXd &eigen_vector)
{
    std::vector<double> v(eigen_vector.data(), eigen_vector.data() + eigen_vector.size());
    return median(v);
}


inline void to_upper(std::string &str)
{
	for(size_t i=0; i<str.size(); i++){
		if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
	}
}


inline int split_string(const std::string& str, std::vector<std::string>& vec_str) {
    vec_str.clear();
    std::istringstream iss(str);
    std::string token;
    while (iss >> token) {
        vec_str.push_back(token);
    }
    return vec_str.size();
}


inline int getMemPeakKB() {
    std::ifstream file("/proc/self/status");
    if (!file.is_open()) return -1;

    std::string line;
    while (std::getline(file, line)) {
        if (line.rfind("VmHWM:", 0) == 0) { // starts with "VmHWM:"
            auto pos = line.find_first_of("0123456789");
            if (pos == std::string::npos) return -1;

            std::istringstream iss(line.substr(pos));
            int value = 0;
            iss >> value;
            return value;
        }
    }
    return -1;
}


inline int getVMPeakKB() {
    std::ifstream file("/proc/self/status");
    if (!file.is_open()) return -1;

    std::string line;
    while (std::getline(file, line)) {
        if (line.rfind("VmPeak:", 0) == 0) { // starts with "VmPeak:"
            auto pos = line.find_first_of("0123456789");
            if (pos == std::string::npos) return -1;

            std::istringstream iss(line.substr(pos));
            int value = 0;
            iss >> value;
            return value;
        }
    }
    return -1;
}
