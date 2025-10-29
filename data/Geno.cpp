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


void Geno::resize(size_t geno_num_, size_t indi_num_) {
    geno_num = geno_num_;
    indi_num = indi_num_;

    words_per_snp = (indi_num + 63) / 64;
    full_bytes = indi_num / 4;
    remainder = indi_num % 4;

    X_A1.resize(geno_num * words_per_snp);
    X_A2.resize(geno_num * words_per_snp);
    X_avg.resize(geno_num);
    X_center_square.resize(geno_num);
    missingness.resize(geno_num);
}


bool Geno::decode_single_genotype(const std::vector<char> &buffer, size_t ref_index, bool swap) {

    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
    const auto& L1 = swap ? LUT_A1_swap : LUT_A1;
    const auto& L2 = swap ? LUT_A2_swap : LUT_A2;
    const size_t base_word_index = ref_index * words_per_snp;

    size_t SNP_sum = 0, SNP_square_sum = 0, not_NA_indi_num = 0, word_index = 0;
    uint64_t wordA1 = 0ULL, wordA2 = 0ULL, bit_index = 0ULL;

    for (size_t j = 0; j < full_bytes; j++) {
        uint8_t b = buffer[j];

        wordA1 |= L1[b] << (bit_index & 63ULL);
        wordA2 |= L2[b] << (bit_index & 63ULL);
        bit_index += 4;

        if ((bit_index & 63ULL) == 0) {
            const auto het = popcnt64(~wordA1 & wordA2); // A1=0, A2=1
            const auto hom = popcnt64(wordA1 & wordA2);  // A1=1, A2=1

            SNP_sum += het + hom * 2;
            SNP_square_sum += het + hom * 4;
            not_NA_indi_num += popcnt64(~wordA1 | wordA2); 

            X_A1[base_word_index + word_index] = wordA1;
            X_A2[base_word_index + word_index] = wordA2;
            wordA1 = wordA2 = 0ULL;
            word_index++;
        }
    }

    // handle remaining genotypes in the last byte
    if (remainder > 0) {
        unsigned char b = buffer[full_bytes];
        const uint64_t keep = (1ULL << remainder) - 1ULL;  // keep only valid genotypes

        wordA1 |= (L1[b] & keep) << (bit_index & 63ULL);
        wordA2 |= (L2[b] & keep) << (bit_index & 63ULL);
        bit_index += remainder;
    }

    // handle leftover partial word, where we may include more than we want
    if ((bit_index & 63ULL) != 0) {
        const uint64_t tail_mask = ((1ULL << (bit_index & 63ULL)) - 1ULL);
        const auto het = popcnt64(~wordA1 & wordA2 & tail_mask); // A1=0, A2=1
        const auto hom = popcnt64(wordA1 & wordA2 & tail_mask);  // A1=1, A2=1

        SNP_sum += het + hom * 2; 
        SNP_square_sum += het + hom * 4;
        not_NA_indi_num += popcnt64((~wordA1 | wordA2) & tail_mask); 

        X_A1[base_word_index + word_index] = wordA1;
        X_A2[base_word_index + word_index] = wordA2;
    }

    if (not_NA_indi_num == 0) return false;

    double temp_not_NA = double(not_NA_indi_num);
    double center_square = SNP_square_sum - SNP_sum * SNP_sum / temp_not_NA;
    if (center_square < 0.5) return false;

    X_avg[ref_index] = SNP_sum / temp_not_NA;
    X_center_square[ref_index] = center_square;
    missingness[ref_index] = 1.0 - temp_not_NA / indi_num;
    return true;
}


double Geno::calc_inner_product(size_t idx1, size_t idx2, bool remove_NA) const
{
    const size_t base_idx1 = idx1 * words_per_snp;
    const size_t base_idx2 = idx2 * words_per_snp;
    uint64_t N_both = 0, S12 = 0, S1 = 0, S2 = 0, S11 = 0, S22 = 0;

    for (size_t w = 0; w < words_per_snp; w++) {
        uint64_t x1_A1 = X_A1[base_idx1 + w], x1_A2 = X_A2[base_idx1 + w];
        uint64_t x2_A1 = X_A1[base_idx2 + w], x2_A2 = X_A2[base_idx2 + w];
        uint64_t both = (~x1_A1 | x1_A2) & (~x2_A1 | x2_A2);

        if (indi_num - w * 64 < 64)
            both &= ((1ULL << (indi_num - w * 64)) - 1ULL);

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
            return (S12 - S1 * S2 / double(N_both)) / std::sqrt(center_square_1) / std::sqrt(center_square_2);
    }
    else
        return (S12 - X_avg[idx1] * S2 - X_avg[idx2] * S1 + X_avg[idx1] * X_avg[idx2] * N_both) / 
            std::sqrt(X_center_square[idx1]) / std::sqrt(X_center_square[idx2]);
}
