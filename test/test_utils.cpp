#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <Eigen/Dense>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/utils.hpp"

TEST_CASE("append_row appends to empty and non-empty matrices") {
    Eigen::MatrixXd m(0, 0);
    Eigen::RowVectorXd r1(3);
    r1 << 1.0, 2.0, 3.0;
    append_row(m, r1);
    CHECK(m.rows() == 1);
    CHECK(m.cols() == 3);
    CHECK(m.row(0).isApprox(r1));

    Eigen::RowVectorXd r2(3);
    r2 << 4.0, 5.0, 6.0;
    append_row(m, r2);
    CHECK(m.rows() == 2);
    CHECK(m.row(1).isApprox(r2));
}

TEST_CASE("append_row validates shape") {
    Eigen::MatrixXd m(0, 0);
    Eigen::VectorXd col(3);
    col << 1.0, 2.0, 3.0;
    CHECK_THROWS_AS(append_row(m, col), std::invalid_argument);

    Eigen::RowVectorXd r1(2);
    r1 << 1.0, 2.0;
    append_row(m, r1);
    Eigen::RowVectorXd r_bad(3);
    r_bad << 1.0, 2.0, 3.0;
    CHECK_THROWS_AS(append_row(m, r_bad), std::invalid_argument);
}

TEST_CASE("append_column appends to empty and non-empty matrices") {
    Eigen::MatrixXd m(0, 0);
    Eigen::VectorXd c1(2);
    c1 << 1.0, 2.0;
    append_column(m, c1);
    CHECK(m.rows() == 2);
    CHECK(m.cols() == 1);
    CHECK(m.col(0).isApprox(c1));

    Eigen::VectorXd c2(2);
    c2 << 3.0, 4.0;
    append_column(m, c2);
    CHECK(m.cols() == 2);
    CHECK(m.col(1).isApprox(c2));
}

TEST_CASE("append_column validates shape") {
    Eigen::MatrixXd m(0, 0);
    Eigen::RowVectorXd row(2);
    row << 1.0, 2.0;
    CHECK_THROWS_AS(append_column(m, row), std::invalid_argument);

    Eigen::VectorXd c1(2);
    c1 << 1.0, 2.0;
    append_column(m, c1);
    Eigen::VectorXd c_bad(3);
    c_bad << 1.0, 2.0, 3.0;
    CHECK_THROWS_AS(append_column(m, c_bad), std::invalid_argument);
}

TEST_CASE("remove_row removes correct index and handles last row") {
    Eigen::MatrixXd m(3, 2);
    m << 1, 2,
         3, 4,
         5, 6;
    remove_row(m, 1);
    CHECK(m.rows() == 2);
    CHECK(m(0, 0) == doctest::Approx(1));
    CHECK(m(1, 0) == doctest::Approx(5));

    remove_row(m, -1);
    CHECK(m.rows() == 1);
    CHECK(m(0, 0) == doctest::Approx(1));
}

TEST_CASE("remove_row on empty matrix is a no-op") {
    Eigen::MatrixXd m(0, 0);
    remove_row(m, -1);
    CHECK(m.rows() == 0);
    CHECK(m.cols() == 0);

    remove_row(m, 0);
    CHECK(m.rows() == 0);
    CHECK(m.cols() == 0);
}

TEST_CASE("remove_column removes correct index and handles last column") {
    Eigen::MatrixXd m(2, 3);
    m << 1, 2, 3,
         4, 5, 6;
    remove_column(m, 1);
    CHECK(m.cols() == 2);
    CHECK(m(0, 0) == doctest::Approx(1));
    CHECK(m(0, 1) == doctest::Approx(3));

    remove_column(m, -1);
    CHECK(m.cols() == 1);
    CHECK(m(0, 0) == doctest::Approx(1));
}

TEST_CASE("remove_column on empty matrix is a no-op") {
    Eigen::MatrixXd m(0, 0);
    remove_column(m, -1);
    CHECK(m.rows() == 0);
    CHECK(m.cols() == 0);

    remove_column(m, 0);
    CHECK(m.rows() == 0);
    CHECK(m.cols() == 0);
}

TEST_CASE("median computes odd/even medians") {
    std::vector<double> v1{3.0, 1.0, 2.0};
    CHECK(median(v1) == doctest::Approx(2.0));

    std::vector<double> v2{4.0, 1.0, 3.0, 2.0};
    CHECK(median(v2) == doctest::Approx(2.5));
}

TEST_CASE("median handles Eigen arrays") {
    Eigen::ArrayXd a(5);
    a << 5.0, 1.0, 3.0, 2.0, 4.0;
    CHECK(median(a) == doctest::Approx(3.0));
}

TEST_CASE("skip_delim and skip_token move the pointer") {
    const char* p = "   token1 token2";
    skip_delim(p);
    CHECK(std::string(p, 6) == "token1");
    skip_token(p);
    skip_delim(p);
    CHECK(std::string(p, 6) == "token2");
}

TEST_CASE("parse_string extracts token and optionally uppercases") {
    const char* p = "  Abc Def";
    std::string out;
    parse_string(p, out, true);
    CHECK(out == "ABC");
    parse_string(p, out, false);
    CHECK(out == "Def");
}

TEST_CASE("parse_num parses numeric tokens") {
    const char* p1 = "  12.5 next";
    double v = 0.0;
    CHECK(parse_num(p1, v) == true);
    CHECK(v == doctest::Approx(12.5));

    const char* p2 = "  12x";
    double v2 = 0.0;
    CHECK(parse_num(p2, v2) == false);
}

TEST_CASE("fast_lookup finds keys and skips invalid entries") {
    std::vector<std::pair<std::string, int>> table = {
        {"a", 0},
        {"b", -1},
        {"c", 2},
    };
    auto it_a = fast_lookup(table, "a");
    CHECK(it_a != table.end());
    CHECK(it_a->second == 0);

    auto it_b = fast_lookup(table, "b");
    CHECK(it_b == table.end());

    auto it_missing = fast_lookup(table, "d");
    CHECK(it_missing == table.end());
}
