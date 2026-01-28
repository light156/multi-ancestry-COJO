#include "include/geno.h"


void Geno::build_LUT() {
    for (int b = 0; b < 256; b++) {
        uint64_t mask1 = 0ULL, mask2 = 0ULL;
        uint64_t mask1_swap = 0ULL, mask2_swap = 0ULL;

        for (int k = 0; k < 4; k++) {
            int code = (b >> (2*k)) & 0x3;

            // normal coding
            switch(code) {
                case 0: // hom major (00) bit1=1, bit2=1
                    mask1 |= 1ULL << k;
                    mask2 |= 1ULL << k;
                    // swapped: hom minor
                    break;
                case 2: // hetero (10) bit1=0, bit2=1
                    mask2 |= 1ULL << k;
                    mask2_swap |= 1ULL << k;
                    break;
                case 3: // hom minor (11) bit1=0, bit2=0
                    // swapped: hom major
                    mask1_swap |= 1ULL << k;
                    mask2_swap |= 1ULL << k;
                    break;
                case 1: // missing (01) bit1=1, bit2=0
                    mask1 |= 1ULL << k;
                    mask1_swap |= 1ULL << k;
                    break;
            }
        }
        
        LUT_bit1[b] = mask1;
        LUT_bit2[b] = mask2;
        LUT_bit1_swap[b] = mask1_swap;
        LUT_bit2_swap[b] = mask2_swap;
    }
}


void Geno::decode_single_genotype(const char* buffer, size_t buffer_size, size_t ref_index, bool swap) {
    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
    const auto& L1 = swap ? LUT_bit1_swap : LUT_bit1;
    const auto& L2 = swap ? LUT_bit2_swap : LUT_bit2;

    auto& SNP_bit1 = X_bit1[ref_index];
    auto& SNP_bit2 = X_bit2[ref_index];
    SNP_bit1.resize(words_per_snp, 0ULL);
    SNP_bit2.resize(words_per_snp, 0ULL);

    int SNP_sum = 0, SNP_square_sum = 0, not_NA_indi_num = 0, byte_index = 0;
    uint64_t word1 = 0ULL, word2 = 0ULL;

    for (size_t w = 0; w < words_per_snp; w++) {
        uint64_t mask = X_mask[w];
        if (mask == 0ULL) {
            byte_index += 16;
            continue;
        }
        
        int bit_index = 0;
        while(bit_index < 64 && byte_index < static_cast<int>(buffer_size)) {    
            uint8_t b = static_cast<uint8_t>(buffer[byte_index]);

            word1 |= (L1[b] << bit_index);
            word2 |= (L2[b] << bit_index);
            bit_index += 4;
            byte_index++;
        }

        const auto het = popcnt64(~word1 & word2 & mask); // bit1=0, bit2=1
        const auto hom = popcnt64(word1 & word2 & mask);  // bit1=1, bit2=1

        SNP_sum += het + hom * 2;
        SNP_square_sum += het + hom * 4;
        not_NA_indi_num += popcnt64((~word1 | word2) & mask);

        SNP_bit1[w] = word1;
        SNP_bit2[w] = word2;
        word1 = word2 = 0ULL;
    }

    if (not_NA_indi_num == 0) return;

    double center_square = SNP_square_sum - double(SNP_sum) * SNP_sum / not_NA_indi_num;
    if (center_square < 0.5) return;

    X_non_NA_indi_num[ref_index] = not_NA_indi_num;
    X_avg[ref_index] = double(SNP_sum) / not_NA_indi_num;
    X_std[ref_index] = std::sqrt(center_square);
}


double Geno::calc_inner_product(size_t idx1, size_t idx2, bool remove_NA) const
{   
    int N_both = 0, S12 = 0, S1 = 0, S2 = 0, S11 = 0, S22 = 0;

    for (size_t w = 0; w < words_per_snp; w++) {
        if (X_mask[w] == 0ULL) continue;

        uint64_t x11 = X_bit1[idx1][w], x12 = X_bit2[idx1][w];
        uint64_t x21 = X_bit1[idx2][w], x22 = X_bit2[idx2][w];
        
        // valid individuals with non-missing genotypes in both SNPs
        uint64_t both = (~x11 | x12) & (~x21 | x22) & X_mask[w]; 
        x11 &= both;
        x12 &= both;
        x21 &= both;
        x22 &= both;

        N_both += popcnt64(both);
        S12 += popcnt64(x11 & x21) + popcnt64(x11 & x22) + popcnt64(x12 & x21) + popcnt64(x12 & x22);

        int n1_1 = popcnt64(~x11 & x12), n1_2 = popcnt64(x11);
        int n2_1 = popcnt64(~x21 & x22), n2_2 = popcnt64(x21);
        
        S1 += n1_1 + n1_2 * 2;
        S2 += n2_1 + n2_2 * 2;

        if (remove_NA) {
            // bit1=0, bit2=1 contributes 1, bit1=1, bit2=1 contributes 4
            S11 += n1_1 + n1_2 * 4;
            S22 += n2_1 + n2_2 * 4;
        }
    }
    
    if (remove_NA) {
        if (N_both == 0) return 0;

        double center_square_1 = S11 - double(S1) * S1 / N_both;
        double center_square_2 = S22 - double(S2) * S2 / N_both;
        if (center_square_1 < 0.5 || center_square_2 < 0.5) return 0;
        
        return (S12 - double(S1) * S2 / N_both) / std::sqrt(center_square_1) / std::sqrt(center_square_2);
    }
    else
        return (S12 - X_avg[idx1] * S2 - X_avg[idx2] * S1 + X_avg[idx1] * X_avg[idx2] * N_both) / X_std[idx1] / X_std[idx2];
}


void Geno::calc_single_genotype_prs(const char* buffer, size_t buffer_size, bool swap, double score_b,
        vector<double>& geno_vec, vector<int>& geno_count_vec, vector<int>& geno_count2_vec)
{
    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
    const auto& L1 = swap ? LUT_bit1_swap : LUT_bit1;
    const auto& L2 = swap ? LUT_bit2_swap : LUT_bit2;

    static thread_local vector<uint64_t> SNP_bit1_tls;
    static thread_local vector<uint64_t> SNP_bit2_tls;
    if (SNP_bit1_tls.size() != words_per_snp) {
        SNP_bit1_tls.assign(words_per_snp, 0ULL);
        SNP_bit2_tls.assign(words_per_snp, 0ULL);
    }

    int SNP_sum = 0, not_NA_indi_num = 0, byte_index = 0;
    uint64_t word1 = 0ULL, word2 = 0ULL;

    for (size_t w = 0; w < words_per_snp; w++) {
        uint64_t mask = X_mask[w];
        if (mask == 0ULL) {
            byte_index += 16;
            continue;
        }
        
        int bit_index = 0;
        while(bit_index < 64 && byte_index < static_cast<int>(buffer_size)) {    
            uint8_t b = static_cast<uint8_t>(buffer[byte_index]);

            word1 |= (L1[b] << bit_index);
            word2 |= (L2[b] << bit_index);
            bit_index += 4;
            byte_index++;
        }

        uint64_t valid = (~word1 | word2) & mask;                
        uint64_t g1 = (~word1 & word2 & mask);   // heterozygous
        uint64_t g2 = (word1 & word2 & mask);   // homozygous alt

        SNP_sum += popcnt64(g1) + popcnt64(g2) * 2;
        not_NA_indi_num += popcnt64(valid);

        SNP_bit1_tls[w] = word1;
        SNP_bit2_tls[w] = word2;

        while (g1) {
            int bit = __builtin_ctzll(g1);
            int n = static_cast<int>(w * 64 + bit);

            geno_vec[n] += score_b;
            geno_count2_vec[n] += 1;

            g1 &= (g1 - 1ULL);
        }

        while (g2) {
            int bit = __builtin_ctzll(g2);
            int n = static_cast<int>(w * 64 + bit);

            geno_vec[n] += score_b * 2.0;
            geno_count2_vec[n] += 2;

            g2 &= (g2 - 1ULL);
        }

        word1 = word2 = 0ULL;
    }

    if (not_NA_indi_num == 0) return;
    double mu = double(SNP_sum) / not_NA_indi_num;

    for (int w = 0; w < words_per_snp; w++) {
        uint64_t miss = (SNP_bit1_tls[w] & ~SNP_bit2_tls[w]) & X_mask[w];

        while (miss) {
            int bit = __builtin_ctzll(miss);
            int n = static_cast<int>(w * 64 + bit);

            geno_vec[n] += mu * score_b;
            geno_count_vec[n] -= 2;

            miss &= (miss - 1ULL);
        }
    }
}