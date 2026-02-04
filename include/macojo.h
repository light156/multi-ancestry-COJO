#pragma once

#include <external/Eigen/Dense>
#include "external/CLI11.hpp"
#include "include/geno.h"
#include "include/config.h"
#include "include/omp_compat.hpp"
#include "include/logger.hpp"
#include "include/utils.hpp"
#include <chrono>

using namespace Eigen;
using namespace std;
using chrono::steady_clock;
using chrono::duration;


// Per-cohort data and computations. Owns the inputs and intermediate matrices
// needed to compute conditional/joint effects under different modes.
class Cohort 
{    
public:
    Cohort(SharedData& s, int cohort_index) : shared(s), cohort_index(cohort_index) {};
    
    // Input readers for this cohort (sumstat/PLINK/LD).
    void read_sumstat();
    void read_frq();
    void read_fam();
    void read_bim();
    void read_bed();
    void read_PLINK_LD();
    // Hidden PRS entry used when --score is provided on the CLI.
    void calc_polygenic_score(int argc, char** argv);

    // Incremental updates to R^{-1} during stepwise selection.
    int calc_R_inv_forward(int append_index);
    int calc_R_inv_backward(int remove_index);

    // Append LD/geno correlations for a new candidate SNP.
    void append_r(const vector<char>& active_mask, int append_index, string mode);
    // Conditional and joint effect estimates for a candidate set.
    void calc_cond_effects(const vector<int>& candidate_SNP, string mode);
    int calc_joint_effects(const vector<int>& candidate_SNP, string mode);
    // Build R^{-1} from scratch for a candidate SNP list.
    int calc_R_inv_from_SNP_list(const vector<int>& SNP_list, string mode); 

// necessary information during calculation
public:
    // sumstat: col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D
    ArrayXXd sumstat;
    MatrixXd r, R_inv_post, R_inv_pre, R_post;
    MatrixXd r_gcta, R_inv_post_gcta, R_inv_pre_gcta; // only used when if_gcta_COJO == true

    ArrayXd beta, beta_var;
    double Vp, R2, previous_R2 = 0.0;
    int cohort_index;

private:
    HyperParams& params = get_params();
    SharedData& shared;

    Geno genotype;
    LDPacked LD_matrix;
    vector<bool> bed_swap_array;
    vector<string> bim_SNP_array, fam_ID_array;
    int valid_indi_num, fam_indi_num;

public:
    // backup for temporary model during backward selection
    struct BackupState {
        MatrixXd r, r_gcta;
        MatrixXd R_inv_pre, R_inv_pre_gcta;
        double previous_R2;
    } backup;

    void save_state() {
        backup.r                 = r;
        backup.r_gcta            = r_gcta;
        backup.R_inv_pre         = R_inv_pre;
        backup.R_inv_pre_gcta    = R_inv_pre_gcta;
        backup.previous_R2       = previous_R2;
    };

    void restore_state() {
        r                 = backup.r;
        r_gcta            = backup.r_gcta;
        R_inv_pre         = backup.R_inv_pre;
        R_inv_pre_gcta    = backup.R_inv_pre_gcta;
        previous_R2       = backup.previous_R2;
    };

    void save_temp_model() {
        R_inv_pre = R_inv_post;
        R_inv_pre_gcta = R_inv_post_gcta;
        previous_R2 = R2;
    };
};


// Orchestrates the multi-cohort analysis workflow and outputs.
class MACOJO
{
public:
    MACOJO() {
        for (int i = 0; i < params.bfile_list.size(); i++) 
            cohorts.emplace_back(shared, i);
    };

    // Read/validate inputs and build common SNP list.
    bool read_input_files();
    // Entry point for joint/conditional/stepwise selection.
    void entry_function();
    
    // Initialize candidate SNP indices from user-provided names.
    void initialize_candidate_SNP(const vector<string>& given_list);
    // Remove collinear SNPs and refresh R^{-1}.
    bool check_candidate_SNP_collinearity(string mode);
    // Stepwise selection loop (forward add, optional backward selection).
    void slct_loop();
    // Inverse-variance meta-analysis across current cohorts.
    void inverse_var_meta(ArrayXd& bma, ArrayXd& se2ma, ArrayXd& abs_zma);

    // Output helpers for results.
    void output_cma(string savename);
    void output_jma(string savename);
    void output_inverse_var_meta(string savename, const vector<pair<int, int>>& SNP_ref_order_pair, bool if_joint);
    void output_ld_matrix(string savename, const vector<pair<int, int>>& SNP_ref_order_pair);
    void output_bad_SNP(string savename);

private:
    HyperParams& params = get_params();
    SharedData shared;

    vector<Cohort> cohorts;
    vector<int> current_list;
    
    vector<int> bad_SNP, candidate_SNP, collinear_SNP, backward_SNP, candidate_SNP_backup, backward_SNP_backup;
    vector<char> active_mask;

    int fixed_candidate_SNP_num = 0;
    
    ArrayXd bC, se2C, abs_zC;
    ArrayXd bJ, se2J, abs_zJ;
};


void skim_fam(string filename, vector<string>& str_list);
void skim_bim(string filename, int chr, vector<string>& str_list);
void skim_SNP(const vector<string>& options, vector<string>& SNP_list);
void set_read_process_output_options(int argc, char** argv);
