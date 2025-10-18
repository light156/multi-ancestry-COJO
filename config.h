#pragma once

#include <map>
#include <string>
#include <vector>


struct HyperParams {
    double threshold = 5e-8;
    double colinear_threshold = 0.9;
    double freq_threshold = 0.01;
    double freq_diff_threshold = 0.2;
    bool if_freq_mode_and = false;
    double R2_incremental_threshold = -1;
    double R2_incremental_threshold_backwards = -1;
    double window_mb = 10;
    int max_iter_num = 10000;

    double iter_colinear_threshold; // 1 / (1 - colinear_threshold)
    int window_size; // in bp, window_mb*1e6

    bool if_gcta_COJO = false;
    bool if_LD_mode = false;
    bool if_cojo_joint = false;
    bool if_MDISA = false;
    bool if_keep_NA = false;
    // bool if_fast_inv = false;
};


struct SharedData {
    uint64_t commonSNP_total_num = 0;
    std::map<std::string,int> commonSNP_index_map;
    std::vector<std::string> A1_ref, A2_ref;
    std::vector<int> SNP_pos_ref, chr_ref;
    std::vector<int> final_commonSNP_index;
    std::vector<std::string> final_commonSNP;
};