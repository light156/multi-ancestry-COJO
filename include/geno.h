#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
using std::vector;
using std::sqrt;
using std::swap;


inline int popcnt64(uint64_t x) {
    return __builtin_popcountll(x);
};


class Geno 
{
public:
    Geno() {
        build_LUT();
    }

    vector<vector<uint64_t>> X_A1, X_A2;
    vector<uint64_t> X_mask;

    vector<double> X_avg, X_center_square;
    vector<int> X_non_NA_indi_num;
    size_t words_per_snp;
    
    void resize(int geno_num) {
        X_A1.resize(geno_num);
        X_A2.resize(geno_num);
        X_avg.resize(geno_num);
        X_center_square.resize(geno_num);
        X_non_NA_indi_num.assign(geno_num, -1);
    };

    void initialize_mask(int valid_indi_num) {
        words_per_snp = (valid_indi_num + 63) / 64;

        X_mask.assign(words_per_snp, ~0ULL);
        if (valid_indi_num % 64 != 0)
            X_mask.back() = ((1ULL << (valid_indi_num % 64)) - 1ULL);
    };

    void decode_single_genotype(const vector<char>& buffer, size_t ref_index, bool swap);
    double calc_inner_product(size_t idx1, size_t idx2, bool remove_NA) const;
    void calc_single_genotype_prs(const vector<char>& buffer, bool swap, double score_b,
        vector<double>& geno_vec, vector<int>& geno_count_vec, vector<int>& geno_count2_vec); 


private:
    uint64_t LUT_A1[256];
    uint64_t LUT_A2[256];
    uint64_t LUT_A1_swap[256];
    uint64_t LUT_A2_swap[256];
    void build_LUT();
};


class LDPacked 
{
public:
    LDPacked() : n_(0) {} // default constructor

    void resize(uint64_t new_n) {
        n_ = new_n;
        data_.assign(new_n * (new_n + 1) / 2, 0.0);
    }

    double& operator()(uint64_t i, uint64_t j) {
        if (i > j) swap(i, j);
        uint64_t idx = i*n_ - i*(i-1)/2 + (j - i);
        return data_[idx];
    }

    double operator()(uint64_t i, uint64_t j) const {
        if (i > j) swap(i, j);
        uint64_t idx = i*n_ - i*(i-1)/2 + (j - i);
        return data_[idx];
    }

    uint64_t size() const { return n_; }

    uint64_t n_;
    vector<double> data_;
};
