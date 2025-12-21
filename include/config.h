#pragma once

#include <string>
#include <vector>
using std::string;
using std::vector;
using std::pair;


struct HyperParams {
    double p = 5e-8;
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


struct SharedData {
    int total_SNP_num = 0;
    vector<pair<string, int>> goodSNP_table;
    vector<string> A1_ref, A2_ref;
    vector<int> SNP_pos_ref;
};
