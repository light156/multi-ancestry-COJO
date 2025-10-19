#pragma once
#include "config.h"
#include "LD.h"
#include "Logger.h"
#include "utils_matrix.h"
#include <omp.h>
#include <map>
#include <vector>
#include <list>
#include <fstream>
#include <Eigen/Dense>
#include <iomanip>
#include <cmath>

using namespace Eigen;
using namespace std;


class Cohort 
{    
public:
    Cohort(const HyperParams& p, SharedData& s) : params(p), shared(s) {};

    void read_sumstat(string cojofile);
    void read_PLINK(string PLINKfile, bool is_ref_cohort);
    void read_PLINK_LD(string PLINKfile, bool is_ref_cohort);
    void skim_fam(string famFile);

    void get_vector_from_bed_matrix(int index, ArrayXd &vec);
    void calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int single_index);

    void calc_conditional_effects();
    bool calc_joint_effects();  
    void calc_adjusted_N();
    bool calc_R_inv();
    void calc_R_inv_gcta(int screened_index);
    bool calc_R_inv_from_SNP_list(const vector<int> &SNP_list, const ArrayXXd &sumstat);
    void save_temp_model();

public:
    // sumstat: col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D, 7:X_avg, 8:X_norm_square
    ArrayXXd sumstat_candidate, sumstat_screened, sumstat_removed;
    
    // X or LD matrix, depend on which mode is used
    vector<uint64_t> X_A1, X_A2;
    LDPacked LD_matrix;
    ArrayXd X_avg, X_square;

    // necessary information during calculation
    MatrixXd r, R_inv_post, R_post, R_inv_pre, R_pre;
    // only used when if_gcta_COJO == true
    MatrixXd r_gcta, R_inv_post_gcta, R_post_gcta, R_inv_pre_gcta, R_pre_gcta;

    // temporary vectors
    ArrayXd X_temp_vec;
    VectorXd r_temp_vec, r_temp_vec_gcta;

    // output 
    ArrayXd conditional_beta, beta, beta_var, output_b, output_se2;
    double Vp, R2, previous_R2 = 0.0;
    int indi_num;

private:
    const HyperParams& params;
    SharedData& shared;

    // backup for temporary model during backward selection
    struct BackupState {
        ArrayXXd sumstat_candidate, sumstat_screened, sumstat_removed;
        MatrixXd r, R_inv_pre, R_pre;
        MatrixXd r_gcta, R_inv_pre_gcta, R_pre_gcta;
        ArrayXd output_b, output_se2;
        double previous_R2;
    };

    BackupState backup;

public:
    void save_state() {
        backup.sumstat_candidate = sumstat_candidate;
        backup.sumstat_screened  = sumstat_screened;
        backup.sumstat_removed   = sumstat_removed;
        backup.r                 = r;
        backup.R_inv_pre         = R_inv_pre;
        backup.R_pre             = R_pre;
        backup.r_gcta            = r_gcta;
        backup.R_inv_pre_gcta    = R_inv_pre_gcta;
        backup.R_pre_gcta        = R_pre_gcta;
        backup.output_b          = output_b;
        backup.output_se2        = output_se2;
        backup.previous_R2       = previous_R2;
    };

    void restore_state() {
        sumstat_candidate = backup.sumstat_candidate;
        sumstat_screened  = backup.sumstat_screened;
        sumstat_removed   = backup.sumstat_removed;
        r                 = backup.r;
        R_inv_pre         = backup.R_inv_pre;
        R_pre             = backup.R_pre;
        r_gcta            = backup.r_gcta;
        R_inv_pre_gcta    = backup.R_inv_pre_gcta;
        R_pre_gcta        = backup.R_pre_gcta;
        output_b          = backup.output_b;
        output_se2        = backup.output_se2;
        previous_R2       = backup.previous_R2;
    };
};


class MACOJO
{
public:
    void read_cojo_PLINK_files(char** filenames, int cohort_num, string extract_file="");
    void read_SNP_only(string filename, vector<string> &SNP_list, bool if_sumstat=false, bool if_bim=false);
    void set_reference_from_bim(string PLINKfile);
    void entry_function(string savename);

    void initialize_main_loop();
    void initialize_MDISA();
    void main_loop();
    void inverse_var_meta_init();
    void inverse_var_meta_conditional();
    void inverse_var_meta_joint();
    void accept_SNP_as_candidate(int screened_index);
    void remove_SNP_from_screened(int screened_index);
    void remove_SNP_from_candidate(int candidate_index);

    void output_results_to_file(string filepath);
    void output_user_hyperparameters();

public:
    HyperParams params;
    SharedData shared;

    vector<Cohort> cohorts;
    vector<int> current_calculation_list;
    vector<int> candidate_SNP, screened_SNP, removed_SNP;
    vector<int> candidate_SNP_backup, screened_SNP_backup, removed_SNP_backup;
    ArrayXd abs_zC, bJ, se2J, abs_zJ;
    
    // user-provided SNP lists
    vector<string> fixed_candidate_SNP;
    int fixed_candidate_SNP_num=0;
};
