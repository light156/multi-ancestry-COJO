#include "include/geno.h"


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


void Geno::decode_single_genotype(const vector<char>& buffer, size_t ref_index, bool swap) {
    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
    const auto& L1 = swap ? LUT_A1_swap : LUT_A1;
    const auto& L2 = swap ? LUT_A2_swap : LUT_A2;

    auto& SNP_A1 = X_A1[ref_index];
    auto& SNP_A2 = X_A2[ref_index];
    SNP_A1.resize(words_per_snp, 0ULL);
    SNP_A2.resize(words_per_snp, 0ULL);

    int SNP_sum = 0, SNP_square_sum = 0, not_NA_indi_num = 0, byte_index = 0;
    uint64_t wordA1 = 0ULL, wordA2 = 0ULL;

    for (size_t w = 0; w < words_per_snp; w++) {
        uint64_t mask = X_mask[w];
        if (mask == 0ULL) {
            byte_index += 16;
            continue;
        }
        
        int bit_index = 0;
        while(bit_index < 64 && byte_index < buffer.size()) {    
            uint8_t b = buffer[byte_index];

            wordA1 |= (L1[b] << bit_index);
            wordA2 |= (L2[b] << bit_index);
            bit_index += 4;
            byte_index++;
        }

        const auto het = popcnt64(~wordA1 & wordA2 & mask); // A1=0, A2=1
        const auto hom = popcnt64(wordA1 & wordA2 & mask);  // A1=1, A2=1

        SNP_sum += het + hom * 2;
        SNP_square_sum += het + hom * 4;
        not_NA_indi_num += popcnt64((~wordA1 | wordA2) & mask);

        SNP_A1[w] = wordA1;
        SNP_A2[w] = wordA2;
        wordA1 = wordA2 = 0ULL;
    }

    if (not_NA_indi_num == 0) return;

    double center_square = SNP_square_sum - double(SNP_sum) * SNP_sum / not_NA_indi_num;
    if (center_square < 0.5) return;

    X_non_NA_indi_num[ref_index] = not_NA_indi_num;
    X_avg[ref_index] = double(SNP_sum) / not_NA_indi_num;
    X_center_square[ref_index] = center_square;
}


double Geno::calc_inner_product(size_t idx1, size_t idx2, bool remove_NA) const
{   
    int N_both = 0, S12 = 0, S1 = 0, S2 = 0, S11 = 0, S22 = 0;

    for (size_t w = 0; w < words_per_snp; w++) {
        if (X_mask[w] == 0ULL) continue;

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
        if (N_both == 0) return 0;

        double center_square_1 = S11 - double(S1) * S1 / N_both;
        double center_square_2 = S22 - double(S2) * S2 / N_both;
        if (center_square_1 < 0.5 || center_square_2 < 0.5) return 0;
        
        return (S12 - double(S1) * S2 / N_both) / sqrt(center_square_1) / sqrt(center_square_2);
    }
    else
        return (S12 - X_avg[idx1] * S2 - X_avg[idx2] * S1 + X_avg[idx1] * X_avg[idx2] * N_both) /
            sqrt(X_center_square[idx1]) / sqrt(X_center_square[idx2]);
}


void Geno::calc_single_genotype_prs(const vector<char>& buffer, bool swap, double score_b,
        vector<double>& geno_vec, vector<int>& geno_count_vec, vector<int>& geno_count2_vec)
{
    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
    const auto& L1 = swap ? LUT_A1_swap : LUT_A1;
    const auto& L2 = swap ? LUT_A2_swap : LUT_A2;

    vector<uint64_t> SNP_A1(words_per_snp, 0ULL);
    vector<uint64_t> SNP_A2(words_per_snp, 0ULL);

    int SNP_sum = 0, not_NA_indi_num = 0, byte_index = 0;
    uint64_t wordA1 = 0ULL, wordA2 = 0ULL;

    for (size_t w = 0; w < words_per_snp; w++) {
        uint64_t mask = X_mask[w];
        if (mask == 0ULL) {
            byte_index += 16;
            continue;
        }
        
        int bit_index = 0;
        while(bit_index < 64 && byte_index < buffer.size()) {    
            uint8_t b = buffer[byte_index];

            wordA1 |= (L1[b] << bit_index);
            wordA2 |= (L2[b] << bit_index);
            bit_index += 4;
            byte_index++;
        }

        uint64_t valid = (~wordA1 | wordA2) & mask;                
        uint64_t g1 = (~wordA1 & wordA2 & mask);   // heterozygous
        uint64_t g2 = (wordA1 & wordA2 & mask);   // homozygous alt

        SNP_sum += popcnt64(g1) + popcnt64(g2) * 2;
        not_NA_indi_num += popcnt64(valid);

        SNP_A1[w] = wordA1;
        SNP_A2[w] = wordA2;

        while (g1) {
            int bit = __builtin_ctzll(g1);
            int n = static_cast<int>(w * 64 + bit);

            geno_vec[n] += score_b;
            geno_count2_vec[n] += 1;

            g1 &= (g1 - 1ULL);
        }

        // g = 2
        while (g2) {
            int bit = __builtin_ctzll(g2);
            int n = static_cast<int>(w * 64 + bit);

            geno_vec[n] += score_b * 2.0;
            geno_count2_vec[n] += 2;

            g2 &= (g2 - 1ULL);
        }

        wordA1 = wordA2 = 0ULL;
    }

    if (not_NA_indi_num == 0) return;
    double mu = double(SNP_sum) / not_NA_indi_num;

    for (int w = 0; w < words_per_snp; w++) {
        uint64_t miss = (SNP_A1[w] & ~SNP_A2[w]) & X_mask[w];

        while (miss) {
            int bit = __builtin_ctzll(miss);
            int n = static_cast<int>(w * 64 + bit);

            geno_vec[n] += mu * score_b;
            geno_count_vec[n] -= 2;

            miss &= (miss - 1ULL);
        }
    }
}