#pragma once

#include "CLI11.hpp"
#include "Logger.h"
#include "Geno.h"
#include "LD.hpp"
#include "utils.hpp"
#include "config.h"
#include "omp_compat.h"
#include <map>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <iomanip>
#include <cmath>

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
    void calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int single_index, VectorXd &r_temp_vec);
    void update_r(const vector<int>& screened_SNP, int append_index);

    // convenient for both iterative selection and effect size calculation
    bool calc_R_inv_from_SNP_list(const vector<int> &SNP_list, bool if_gcta_COJO, bool if_remove_NA);
    bool calc_joint_effects(bool if_gcta_COJO);

// necessary information during calculation
public:
    // sumstat: col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D
    ArrayXXd sumstat, sumstat_candidate;

    MatrixXd r, R_inv_post, R_inv_pre, R_post;
    MatrixXd r_gcta, R_inv_post_gcta, R_inv_pre_gcta; // only used when if_gcta_COJO == true

    ArrayXd beta, beta_var;
    double Vp, R2, previous_R2;
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

    // backup for temporary model during backward selection
    struct BackupState {
        ArrayXXd sumstat_candidate;
        MatrixXd r, R_inv_pre;
        MatrixXd r_gcta, R_inv_pre_gcta;
        double previous_R2;
    };

    BackupState backup;

public:
    void save_state() {
        backup.sumstat_candidate = sumstat_candidate;
        backup.r                 = r;
        backup.R_inv_pre         = R_inv_pre;
        backup.r_gcta            = r_gcta;
        backup.R_inv_pre_gcta    = R_inv_pre_gcta;
        backup.previous_R2       = previous_R2;
    };

    void restore_state() {
        sumstat_candidate = backup.sumstat_candidate;
        r                 = backup.r;
        R_inv_pre         = backup.R_inv_pre;
        r_gcta            = backup.r_gcta;
        R_inv_pre_gcta    = backup.R_inv_pre_gcta;
        previous_R2       = backup.previous_R2;
        backup            = BackupState(); // clear backup after restore
    };

    void save_temp_model() 
    {   
        if (sumstat_candidate.rows() == 1) {
            R_inv_pre = MatrixXd::Identity(1,1);
            R_inv_pre_gcta = MatrixXd::Identity(1,1) / sumstat_candidate(0,6); // only used when if_gcta_COJO == true
            previous_R2 = -1; // only used for if_gcta_COJO == false;
        } else {
            R_inv_pre = R_inv_post;
            R_inv_pre_gcta = R_inv_post_gcta; // only used when if_gcta_COJO == true
            previous_R2 = R2; // only used for if_gcta_COJO == false
        }
    };
};


class MACOJO
{
public:
    int set_read_process_output_options(int argc, char** argv);
    void read_cojo_PLINK_files();
    void entry_function();
    void output_cojo_cond(string savename);
    void output_cojo_joint(string savename);

    void main_loop();
    void inverse_var_meta_conditional();
    void inverse_var_meta_joint();

private:
    HyperParams params;
    SharedData shared;

    vector<Cohort> cohorts;
    vector<int> current_list;
    vector<int> bad_SNP, screened_SNP;
    vector<int> candidate_SNP, removed_SNP, candidate_SNP_backup, removed_SNP_backup;
    int fixed_candidate_SNP_num = 0;
    
    ArrayXd bC, se2C, abs_zC;
    ArrayXd bJ, se2J, abs_zJ;
};


void skim_file(string filename, vector<string> &SNP_list, bool header, bool first_column, bool second_column);