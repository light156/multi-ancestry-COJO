#pragma once

#include <Eigen/Dense>
#include "CLI11.hpp"

#include "Geno.h"
#include "LD.hpp"
#include "config.h"
#include "omp_compat.hpp"
#include "logger.hpp"
#include "utils.hpp"
#include <chrono>

using namespace Eigen;
using namespace std;


class Cohort 
{    
public:
    Cohort(const HyperParams& p, SharedData& s, int cohort_index) : 
        params(p), shared(s), cohort_index(cohort_index) {};

    void read_sumstat();
    void read_frq();
    void read_fam();
    void read_bim();
    void read_bed();
    void read_PLINK_LD();

    bool calc_R_inv_forward(int append_index);
    void calc_R_inv_backward(int remove_index);

    void append_r(const vector<int>& SNP_list, int append_index, string mode);
    void calc_cond_effects(const vector<int>& candidate_SNP, string mode);
    bool calc_joint_effects(const vector<int>& candidate_SNP, string mode);
    int calc_R_inv_from_SNP_list(const vector<int> &SNP_list, string mode); 

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
    const HyperParams& params;
    SharedData& shared;

    Geno genotype;
    LDPacked LD_matrix;
    vector<bool> swap_array;
    vector<string> bim_SNP_list;
    vector<int> fam_keep_list;
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


class MACOJO
{
public:
    int set_read_process_output_options(int argc, char** argv);
    void read_input_files();
    void initialize_candidate_SNP(string filename);
    void slct_loop();
    void inverse_var_meta(ArrayXd &bma, ArrayXd &se2ma, ArrayXd &abs_zma);
    bool check_candidate_SNP_collinearity(string mode);

    void entry_function();
    void output_cma(string savename);
    void output_jma(string savename);
    void output_inverse_var_meta(string savename, char mode, const map<int, int>& SNP_ref_order_pair);
    void output_ld_matrix(string savename, const vector<int>& ordered_candidate, const Cohort& c);

private:
    HyperParams params;
    SharedData shared;

    vector<Cohort> cohorts;
    vector<int> current_list;
    vector<int> bad_SNP, screened_SNP;
    vector<int> candidate_SNP, collinear_SNP, backward_SNP, candidate_SNP_backup, backward_SNP_backup;
    int fixed_candidate_SNP_num = 0;
    
    ArrayXd bC, se2C, abs_zC;
    ArrayXd bJ, se2J, abs_zJ;
};


void skim_file(string filename, vector<string> &SNP_list, bool header, bool first_column, bool second_column);