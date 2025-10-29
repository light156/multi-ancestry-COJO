#pragma once

#include "CLI11.hpp"
#include "Logger.h"
#include "Geno.h"
#include "LD.hpp"
#include "utils_matrix.hpp"
#include "config.h"
#include "omp_compat.h"
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
    void read_PLINK_bed(string PLINKfile, bool is_ref_cohort);
    void read_PLINK_LD(string PLINKfile, bool is_ref_cohort);
    void skim_fam(string PLINKfile);
    void skim_frq(string PLINKfile);

    bool calc_R_inverse_forward(int append_index);
    void calc_R_inverse_backward(int remove_index);
    bool calc_R_inv_from_SNP_list(const vector<int> &SNP_list);
    
    void calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int single_index, VectorXd &r_temp_vec);
    void calc_conditional_effects();
    bool calc_joint_effects();
    void save_temp_model();

public:
    // sumstat: col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D
    ArrayXXd sumstat, sumstat_candidate;
    Geno genotype;
    LDPacked LD_matrix;
    vector<int> bim_index_list; // index mapping to shared common SNP list, -1 if not in common list
    vector<bool> swap_array;

    // necessary information during calculation
    MatrixXd r, R_inv_post, R_inv_pre, R_post;
    MatrixXd r_gcta, R_inv_post_gcta, R_inv_pre_gcta; // only used when if_gcta_COJO == true

    // output 
    ArrayXd conditional_beta, beta, beta_var, output_b, output_se2;
    double Vp, R2, previous_R2;
    int indi_num;

private:
    const HyperParams& params;
    SharedData& shared;

    // backup for temporary model during backward selection
    struct BackupState {
        ArrayXXd sumstat_candidate;
        MatrixXd r, R_inv_pre;
        MatrixXd r_gcta, R_inv_pre_gcta;
        ArrayXd output_b, output_se2;
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
        backup.output_b          = output_b;
        backup.output_se2        = output_se2;
        backup.previous_R2       = previous_R2;
    };

    void restore_state() {
        sumstat_candidate = backup.sumstat_candidate;
        r                 = backup.r;
        R_inv_pre         = backup.R_inv_pre;
        r_gcta            = backup.r_gcta;
        R_inv_pre_gcta    = backup.R_inv_pre_gcta;
        output_b          = backup.output_b;
        output_se2        = backup.output_se2;
        previous_R2       = backup.previous_R2;
        backup            = BackupState(); // clear backup after restore
    };
};


class MACOJO
{
public:
    void read_SNP_only(string filename, vector<string> &SNP_list, bool if_sumstat=false, bool if_bim=false);
    void read_cojo_PLINK_files();
    void set_reference_from_bim();
    void output_results_to_file(string savename);
    void entry_function();
    int set_read_process_output_options(int argc, char** argv);

    void initialize_from_fixed_candidate();
    void main_loop();
    void inverse_var_meta_conditional();
    void inverse_var_meta_joint();
    void accept_SNP_as_candidate(int index);
    void remove_SNP_from_candidate(int candidate_index);

private:
    HyperParams params;
    SharedData shared;
    
    vector<Cohort> cohorts;
    vector<int> current_list;
    vector<int> good_SNP, bad_SNP, screened_SNP;
    vector<int> candidate_SNP, removed_SNP, candidate_SNP_backup, removed_SNP_backup;
    ArrayXd abs_zC, bJ, se2J, abs_zJ;

    vector<string> fixed_candidate_SNP;
    int fixed_candidate_SNP_num = 0;

    // filepaths
    string extract_file="", fixedSNP_file="", output_name="";
    vector<string> bfile_list, cojo_file_list;
};


inline double median(const std::vector<double> &v)
{
    int size = v.size();
	if (size == 0) 
        throw std::invalid_argument("median: empty vector");

    std::vector<double> b(v);
    std::sort(b.begin(), b.end());
	return (size%2==1) ? b[(size-1)/2] : (b[size/2]+b[size/2-1])/2;
}


inline double median(const ArrayXd &eigen_vector)
{
    std::vector<double> v(eigen_vector.data(), eigen_vector.data() + eigen_vector.size());
    return median(v);
}
