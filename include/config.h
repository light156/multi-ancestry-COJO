#pragma once

#include <map>
#include <string>
#include <vector>


struct HyperParams {
    double p = 5e-8;
    double collinear = 0.9;
    double maf = 0.01;
    double missingness = 1;
    double diff_freq = 0.2;
    bool if_freq_mode_and = false;
    double R2_threshold = -1;
    double R2back_threshold = -1;
    double window_kb = 10000;
    int max_iter_num = 10000;

    double iter_collinear_threshold; // 1 / (1 - collinear)
    int window_size; // in bp, window_kb*1e3

    string effect_size_mode; // "GCTA", "removeNA", "imputeNA"
    bool if_gcta_COJO, if_remove_NA;

    bool if_joint_mode = false;
    bool if_cond_mode = false;
    bool if_MDISA = false;
    bool if_LD_mode = false;

    // filepaths
    std::vector<std::string> bfile_list, cojo_file_list, keep_file_list, remove_file_list;
    string output_name, extract_file, exclude_file, fixedSNP_file, cond_file;
};


struct SharedData {
    int total_SNP_num = 0;
    std::map<std::string, int> goodSNP_index_map;
    std::vector<std::string> SNP_ref, A1_ref, A2_ref;
    std::vector<int> SNP_pos_ref, chr_ref;
};
