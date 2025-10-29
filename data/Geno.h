#pragma once
#include <vector>
#include <cstdint>
#include <cmath>


class Geno 
{
public:
    Geno() {
        build_LUT();
    }

    std::vector<uint64_t> X_A1, X_A2;
    std::vector<double> X_avg, X_center_square, missingness;
    
    size_t indi_num, geno_num;
    size_t words_per_snp, full_bytes, remainder;
    
    void resize(size_t geno_num_, size_t indi_num_);
    bool decode_single_genotype(const std::vector<char> &buffer, size_t ref_index, bool swap);
    double calc_inner_product(size_t idx1, size_t idx2, bool remove_NA) const;

private:
    uint64_t LUT_A1[256];
    uint64_t LUT_A2[256];
    uint64_t LUT_A1_swap[256];
    uint64_t LUT_A2_swap[256];
    void build_LUT();
};


inline int popcnt64(uint64_t x) {
    return __builtin_popcountll(x);
};