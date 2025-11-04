#pragma once
#include <vector>
#include <cstdint>
#include <cmath>
using std::vector;
using std::sqrt;


class Geno 
{
public:
    Geno() {
        build_LUT();
        build_LUT_single();
    }

    vector<vector<uint64_t>> X_A1, X_A2;
    vector<double> X_avg, X_center_square;
    vector<uint64_t> X_mask;
    size_t words_per_snp;
    
    void resize(int geno_num) {
        X_A1.resize(geno_num);
        X_A2.resize(geno_num);
        X_avg.resize(geno_num);
        X_center_square.resize(geno_num);
    };

    void initialize_mask(int valid_indi_num) {
        words_per_snp = (valid_indi_num + 63) / 64;

        X_mask.assign(words_per_snp, ~0ULL);
        if (valid_indi_num % 64 != 0)
            X_mask.back() = ((1ULL << (valid_indi_num % 64)) - 1ULL);
    };

    int decode_single_genotype(const vector<char>& buffer, size_t ref_index, bool swap);
    int decode_single_genotype(const vector<char>& buffer, size_t ref_index, bool swap, const vector<int>& keep_list);
    double calc_inner_product(size_t idx1, size_t idx2, bool remove_NA) const;

private:
    uint64_t LUT_A1[256];
    uint64_t LUT_A2[256];
    uint64_t LUT_A1_swap[256];
    uint64_t LUT_A2_swap[256];
    void build_LUT();

    uint64_t LUT_A1_single[256][4];
    uint64_t LUT_A2_single[256][4];
    uint64_t LUT_A1_swap_single[256][4];
    uint64_t LUT_A2_swap_single[256][4];
    void build_LUT_single();
};


inline int popcnt64(uint64_t x) {
    return __builtin_popcountll(x);
};
