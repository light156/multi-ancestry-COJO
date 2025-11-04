#include "Geno.h"


void Geno::build_LUT() {
    for (int b = 0; b < 256; b++) {
        uint64_t A1mask = 0ULL, A2mask = 0ULL;
        uint64_t A1mask_swap = 0ULL, A2mask_swap = 0ULL;

        for (int k = 0; k < 4; k++) {
            int code = (b >> (2*k)) & 0x3;

            // normal coding
            switch(code) {
                case 0: // hom major (00) A1=1, A2=1
                    A1mask |= 1ULL << k;
                    A2mask |= 1ULL << k;
                    // swapped: hom minor
                    break;
                case 2: // hetero (10) A1=0, A2=1
                    A2mask |= 1ULL << k;
                    A2mask_swap |= 1ULL << k;
                    break;
                case 3: // hom minor (11) A1=0, A2=0
                    // swapped: hom major
                    A1mask_swap |= 1ULL << k;
                    A2mask_swap |= 1ULL << k;
                    break;
                case 1: // missing (01) A1=1, A2=0
                    A1mask |= 1ULL << k;
                    A1mask_swap |= 1ULL << k;
                    break;
            }
        }
        
        LUT_A1[b] = A1mask;
        LUT_A2[b] = A2mask;
        LUT_A1_swap[b] = A1mask_swap;
        LUT_A2_swap[b] = A2mask_swap;
    }
}


void Geno::build_LUT_single() {
    for (int b = 0; b < 256; b++) {
        for (int k = 0; k < 4; k++) {
            int code = (b >> (2*k)) & 0x3;

            uint64_t A1 = 0, A2 = 0;
            uint64_t A1s = 0, A2s = 0;

            switch (code) {
                case 0: // hom major (00) A1=1, A2=1
                    A1 = 1; A2 = 1; break;
                case 2: // hetero (10) A1=0, A2=1
                    A2 = 1; A2s = 1; break;
                case 3: // hom minor (11) A1=0, A2=0
                    // swapped: hom major
                    A1s = 1; A2s = 1; break;
                case 1: // missing (01) A1=1, A2=0
                    A1 = 1; A1s = 1; break;
            }

            LUT_A1_single[b][k] = A1;
            LUT_A2_single[b][k] = A2;
            LUT_A1_swap_single[b][k] = A1s;
            LUT_A2_swap_single[b][k] = A2s;
        }
    }
}


int Geno::decode_single_genotype(const vector<char>& buffer, size_t ref_index, bool swap) {

    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
    const auto& L1 = swap ? LUT_A1_swap : LUT_A1;
    const auto& L2 = swap ? LUT_A2_swap : LUT_A2;

    auto& SNP_A1 = X_A1[ref_index];
    auto& SNP_A2 = X_A2[ref_index];
    SNP_A1.resize(words_per_snp);
    SNP_A2.resize(words_per_snp);

    int SNP_sum = 0, SNP_square_sum = 0, not_NA_indi_num = 0, word_index = 0, bit_index = 0;
    uint64_t wordA1 = 0ULL, wordA2 = 0ULL;

    for (int j = 0; j < buffer.size(); j++) {    
        uint8_t b = buffer[j];

        wordA1 |= (L1[b] << bit_index);
        wordA2 |= (L2[b] << bit_index);
        bit_index += 4;

        if (bit_index == 64 || j == buffer.size() - 1) {
            const auto mask = X_mask[word_index]; // mask for valid individuals in this word
            const auto het = popcnt64(~wordA1 & wordA2 & mask); // A1=0, A2=1
            const auto hom = popcnt64(wordA1 & wordA2 & mask);  // A1=1, A2=1

            SNP_sum += het + hom * 2;
            SNP_square_sum += het + hom * 4;
            not_NA_indi_num += popcnt64((~wordA1 | wordA2) & mask);

            SNP_A1[word_index] = wordA1;
            SNP_A2[word_index] = wordA2;
            wordA1 = wordA2 = 0ULL;
            bit_index = 0;
            word_index++;
        }
    }

    if (not_NA_indi_num == 0) return false;

    double center_square = SNP_square_sum - double(SNP_sum) * SNP_sum / not_NA_indi_num;
    if (center_square < 0.5) return false;

    X_avg[ref_index] = double(SNP_sum) / not_NA_indi_num;
    X_center_square[ref_index] = center_square;
    return not_NA_indi_num;
}


int Geno::decode_single_genotype(const vector<char>& buffer, size_t ref_index, bool swap, const vector<int>& keep_list) {

    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
    const auto& L1 = swap ? LUT_A1_swap_single : LUT_A1_single;
    const auto& L2 = swap ? LUT_A2_swap_single : LUT_A2_single;

    auto& SNP_A1 = X_A1[ref_index];
    auto& SNP_A2 = X_A2[ref_index];
    SNP_A1.resize(words_per_snp);
    SNP_A2.resize(words_per_snp);

    int SNP_sum = 0, SNP_square_sum = 0, not_NA_indi_num = 0, word_index = 0, bit_index = 0;
    uint64_t wordA1 = 0ULL, wordA2 = 0ULL;

    for (int j = 0; j < keep_list.size(); j++) {
        uint8_t b = buffer[keep_list[j] >> 2];

        wordA1 |= (L1[b][keep_list[j] & 3] << bit_index);
        wordA2 |= (L2[b][keep_list[j] & 3] << bit_index);
        bit_index++;

        if (bit_index == 64 || j == keep_list.size() - 1) {
            const auto mask = X_mask[word_index]; // mask for valid individuals in this word
            const auto het = popcnt64(~wordA1 & wordA2 & mask); // A1=0, A2=1
            const auto hom = popcnt64(wordA1 & wordA2 & mask);  // A1=1, A2=1

            SNP_sum += het + hom * 2;
            SNP_square_sum += het + hom * 4;
            not_NA_indi_num += popcnt64((~wordA1 | wordA2) & mask);

            SNP_A1[word_index] = wordA1;
            SNP_A2[word_index] = wordA2;
            wordA1 = wordA2 = 0ULL;
            bit_index = 0;
            word_index++;
        }
    }

    if (not_NA_indi_num == 0) return false;

    double center_square = SNP_square_sum - double(SNP_sum) * SNP_sum / not_NA_indi_num;
    if (center_square < 0.5) return false;

    X_avg[ref_index] = double(SNP_sum) / not_NA_indi_num;
    X_center_square[ref_index] = center_square;
    return not_NA_indi_num;
}


double Geno::calc_inner_product(size_t idx1, size_t idx2, bool remove_NA) const
{   
    uint64_t N_both = 0, S12 = 0, S1 = 0, S2 = 0, S11 = 0, S22 = 0;

    for (size_t w = 0; w < words_per_snp; w++) {
        uint64_t x1_A1 = X_A1[idx1][w], x1_A2 = X_A2[idx1][w];
        uint64_t x2_A1 = X_A1[idx2][w], x2_A2 = X_A2[idx2][w];
        
        // valid individuals with non-missing genotypes in both SNPs
        uint64_t both = (~x1_A1 | x1_A2) & (~x2_A1 | x2_A2) & X_mask[w]; 
        x1_A1 &= both;
        x1_A2 &= both;
        x2_A1 &= both;
        x2_A2 &= both;

        N_both += popcnt64(both);
        S1 += popcnt64(x1_A1) + popcnt64(x1_A2);
        S2 += popcnt64(x2_A1) + popcnt64(x2_A2);
        // (x1_A1+x1_A2) * (x2_A1+x2_A2) = x1_A1&x2_A1 + x1_A1&x2_A2 + x1_A2&x2_A1 + x1_A2&x2_A2
        S12 += popcnt64(x1_A1 & x2_A1) + popcnt64(x1_A1 & x2_A2) + popcnt64(x1_A2 & x2_A1) + popcnt64(x1_A2 & x2_A2);

        if (remove_NA) {
            // A1=0, A2=1 contributes 1, A1=1, A2=1 contributes 4
            S11 += popcnt64(~x1_A1 & x1_A2) + popcnt64(x1_A1 & x1_A2) * 4;
            S22 += popcnt64(~x2_A1 & x2_A2) + popcnt64(x2_A1 & x2_A2) * 4;
        }
    }
    
    if (remove_NA) {
        double center_square_1 = S11 - S1 * S1 / double(N_both);
        double center_square_2 = S22 - S2 * S2 / double(N_both);
        if (N_both == 0 || center_square_1 < 0.5 || center_square_2 < 0.5) 
            return 0;
        else
            return (S12 - S1 * S2 / double(N_both)) / sqrt(center_square_1) / sqrt(center_square_2);
    }
    else
        return (S12 - X_avg[idx1] * S2 - X_avg[idx2] * S1 + X_avg[idx1] * X_avg[idx2] * N_both) /
            sqrt(X_center_square[idx1]) / sqrt(X_center_square[idx2]);
}
