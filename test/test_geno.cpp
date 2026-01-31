#include "doctest.h"

#include <cstdint>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "../include/geno.h"

static uint8_t encode_geno_2bit(int g) {
    // Encoding expected by decode_single_genotype (from build_LUT logic):
    // 00 -> hom (value 2), 10 -> hetero (value 1), 11 -> hom (value 0), 01 -> missing
    switch (g) {
        case 2: return 0b00;
        case 1: return 0b10;
        case 0: return 0b11;
        default: return 0b01;
    }
}

static std::vector<char> build_snp_buffer(const std::vector<int>& genotypes) {
    const size_t n = genotypes.size();
    const size_t bytes_per_snp = (n + 3) / 4;
    std::vector<char> buffer(bytes_per_snp, 0);

    for (size_t i = 0; i < n; ++i) {
        const size_t byte_idx = i / 4;
        const size_t shift = (i % 4) * 2;
        const uint8_t code = encode_geno_2bit(genotypes[i]);
        buffer[byte_idx] |= static_cast<char>(code << shift);
    }

    return buffer;
}

static double pearson_corr(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    Eigen::VectorXd xc = x.array() - x.mean();
    Eigen::VectorXd yc = y.array() - y.mean();
    double num = xc.dot(yc);
    double den = std::sqrt(xc.squaredNorm() * yc.squaredNorm());
    if (den == 0.0) return 0.0;
    return num / den;
}

static double expected_inner_remove_na(
    const std::vector<int>& g1,
    const std::vector<int>& g2)
{
    std::vector<double> v1, v2;
    v1.reserve(g1.size());
    v2.reserve(g2.size());
    for (size_t i = 0; i < g1.size(); ++i) {
        if (g1[i] < 0 || g2[i] < 0) continue;
        v1.push_back(static_cast<double>(g1[i]));
        v2.push_back(static_cast<double>(g2[i]));
    }
    if (v1.empty()) return 0.0;
    Eigen::Map<Eigen::VectorXd> x1(v1.data(), static_cast<Eigen::Index>(v1.size()));
    Eigen::Map<Eigen::VectorXd> x2(v2.data(), static_cast<Eigen::Index>(v2.size()));
    if ((x1.array() - x1.mean()).square().sum() < 0.5 ||
        (x2.array() - x2.mean()).square().sum() < 0.5) return 0.0;
    return pearson_corr(x1, x2);
}

static double expected_inner_impute_na(
    const std::vector<int>& g1,
    const std::vector<int>& g2)
{
    std::vector<double> v1_all, v2_all;
    v1_all.reserve(g1.size());
    v2_all.reserve(g2.size());
    for (int v : g1) if (v >= 0) v1_all.push_back(static_cast<double>(v));
    for (int v : g2) if (v >= 0) v2_all.push_back(static_cast<double>(v));
    if (v1_all.empty() || v2_all.empty()) return 0.0;

    Eigen::Map<Eigen::VectorXd> x1_all(v1_all.data(), static_cast<Eigen::Index>(v1_all.size()));
    Eigen::Map<Eigen::VectorXd> x2_all(v2_all.data(), static_cast<Eigen::Index>(v2_all.size()));
    const double cs1 = (x1_all.array() - x1_all.mean()).square().sum();
    const double cs2 = (x2_all.array() - x2_all.mean()).square().sum();
    if (cs1 < 0.5 || cs2 < 0.5) return 0.0;

    Eigen::VectorXd x1_imputed(g1.size());
    Eigen::VectorXd x2_imputed(g2.size());
    for (size_t i = 0; i < g1.size(); ++i) {
        x1_imputed[i] = (g1[i] < 0) ? x1_all.mean() : static_cast<double>(g1[i]);
        x2_imputed[i] = (g2[i] < 0) ? x2_all.mean() : static_cast<double>(g2[i]);
    }
    return pearson_corr(x1_imputed, x2_imputed);
}

TEST_CASE("decode_single_genotype handles random genotypes with missing") {
    const int n_individuals = 100;
    const int n_snps = 5;

    Geno g;
    g.resize(n_snps);
    g.initialize_mask(n_individuals);

    std::mt19937 rng(1337);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::uniform_int_distribution<int> geno_dist(0, 2);

    for (int snp = 0; snp < n_snps; ++snp) {
        std::vector<int> genotypes;
        genotypes.reserve(n_individuals);

        int non_na = 0;
        double sum = 0.0;
        double sumsq = 0.0;

        for (int i = 0; i < n_individuals; ++i) {
            int gval = geno_dist(rng);
            if (u(rng) < 0.1) gval = -1; // ~10% missing
            genotypes.push_back(gval);
            if (gval >= 0) {
                non_na++;
                sum += gval;
                sumsq += gval * gval;
            }
        }

        const std::vector<char> buffer = build_snp_buffer(genotypes);
        g.decode_single_genotype(buffer.data(), buffer.size(), snp, false);

        if (non_na == 0) {
            CHECK(g.X_non_NA_indi_num[snp] == -1);
            continue;
        }

        const double mean = sum / non_na;
        const double center_square = sumsq - sum * sum / non_na;

        if (center_square < 0.5) {
            CHECK(g.X_non_NA_indi_num[snp] == -1);
        } else {
            CHECK(g.X_non_NA_indi_num[snp] == non_na);
            CHECK(g.X_avg[snp] == doctest::Approx(mean));
            CHECK(g.X_std[snp] == doctest::Approx(std::sqrt(center_square)));
        }
    }
}

TEST_CASE("calc_inner_product matches reference formulas (remove_NA true/false)") {
    const int n_individuals = 20;
    Geno g;
    g.resize(2);
    g.initialize_mask(n_individuals);

    std::vector<int> g1{0, 1, 2, -1, 2, 0, 1, 2, -1, 0, 1, 2, 0, 2, 1, -1, 2, 0, 1, 2};
    std::vector<int> g2{2, 1, 0,  2, -1, 1, 0, 2,  1, -1, 2, 0, 1, 0, 2,  1, -1, 2, 0, 1};

    std::vector<char> b1 = build_snp_buffer(g1);
    std::vector<char> b2 = build_snp_buffer(g2);

    g.decode_single_genotype(b1.data(), b1.size(), 0, false);
    g.decode_single_genotype(b2.data(), b2.size(), 1, false);

    const double expected_remove = expected_inner_remove_na(g1, g2);
    const double expected_keep = expected_inner_impute_na(g1, g2);

    const double got_remove = g.calc_inner_product(0, 1, true);
    const double got_keep = g.calc_inner_product(0, 1, false);

    CHECK(got_remove == doctest::Approx(expected_remove));
    CHECK(got_keep == doctest::Approx(expected_keep));

    CHECK(g.calc_inner_product(0, 0, true) == doctest::Approx(1.0));
    CHECK(g.calc_inner_product(1, 1, true) == doctest::Approx(1.0));
    CHECK(g.calc_inner_product(0, 0, false) == doctest::Approx(1.0));
    CHECK(g.calc_inner_product(1, 1, false) == doctest::Approx(1.0));
}
