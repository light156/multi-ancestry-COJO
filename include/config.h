#pragma once

#include <string>
#include <vector>
using std::string;
using std::vector;
using std::pair;


// Global configuration derived from CLI arguments.
struct HyperParams {
    double p_value = 5e-8;
    double collinear = 0.9;
    double maf = 0.01;
    double missingness = 1;
    double diff_freq = 0.2;
    bool if_freq_mode_and = false;
    bool if_infoscore = false;
    double R2_threshold = -1;
    double R2back_threshold = -1;
    double window_kb = 10000;
    int max_iter_num = 10000;
    int thread_num = 1;
    int bed_block_mb = 64;

    double iter_collinear_threshold; // 1 / (1 - collinear)
    double window_size; // in bp, window_kb*1e3

    vector<int> chr_list = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
    int curr_chr = -1;

    string slct_mode, effect_size_mode; // "GCTA", "removeNA", "imputeNA"

    bool if_joint_mode = false;
    bool if_cond_mode = false;
    bool if_MDISA = false;
    bool if_LD_mode = false;
    bool if_output_all = false;

    // input data filepaths
    vector<string> bfile_list, cojo_file_list, keep_file_list, remove_file_list;
    string output_name;

    // SNP lists
    vector<string> extract_options, exclude_options, fix_options, cojo_cond_options;
    vector<string> extract_SNPs, exclude_SNPs, fix_SNPs, cojo_cond_SNPs;
};


inline HyperParams& get_params() {
    static HyperParams p;
    return p;
}


enum class BadSnpReason : uint8_t {
    FreqDiffTooLarge    = 0,
    HighMissingness     = 1,
    AlleleMismatch      = 2,
    BpMismatch          = 3,
    InvalidNumericValue = 4,
    RareVariant         = 5,
};


inline const char* reason_to_string(BadSnpReason r) {
    switch (r) {
        case BadSnpReason::FreqDiffTooLarge:    return "Allele frequency too different between sumstat and genotype data";
        case BadSnpReason::HighMissingness:     return "Genotype missingness too high or all values identical in bedfile";
        case BadSnpReason::AlleleMismatch:      return "A1 and A2 different from ref BIM file";
        case BadSnpReason::BpMismatch:          return "SNP position different from ref BIM file";
        case BadSnpReason::InvalidNumericValue: return "Invalid numeric value in sumstat file";
        case BadSnpReason::RareVariant:         return "Allele frequency below MAF threshold";
    }
    return "Unknown";
}


// Shared SNP references across cohorts (aligned to cohort 1).
struct SharedData {
    vector<pair<string, int>> goodSNP_table;
    vector<string> A1_ref, A2_ref;
    vector<int> SNP_pos_ref;
    vector<pair<BadSnpReason, string>> bad_SNP_dict;
    vector<int> bp_order;
};
