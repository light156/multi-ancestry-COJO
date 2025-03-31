#ifndef _TACOJO_H
#define _TACOJO_H

#include "Logger.h"
#include <omp.h>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <bitset>
#include <Eigen/Dense>
#include <unsupported/Eigen/SpecialFunctions>
#include <iomanip>
#include <numeric>

using namespace Eigen;
using namespace std;


void to_upper(string &str);
int split_string(const string &str, vector<string> &vec_str, string separator=" ,\t;\n");
double median(const ArrayXd &eigen_vector);

void append_row(ArrayXXd &matrix, const ArrayXXd &vector);
void append_row(MatrixXd &matrix, const MatrixXd &vector);
void append_column(MatrixXd &matrix, const MatrixXd &vector);
void remove_row(ArrayXXd &matrix, int index=-1);
void remove_row(MatrixXd &matrix, int index=-1);
void remove_column(MatrixXd &matrix, int index=-1);


class Cohort 
{    
public:
    void skim_cojo(string cojofile);
    void skim_PLINK(string PLINKfile);
    void read_sumstat(string cojofile);
    void read_PLINK(string PLINKfile, bool is_ref_cohort);
    
    void get_vector_from_bed_matrix(int index, VectorXd &vec);
    void calc_conditional_effects();
    bool calc_joint_effects(const ArrayXXd &sumstat_temp, bool flag, double iter_colinear_threshold=0);    
    void calc_Vp(ArrayXXd &sumstat_temp);
    void calc_R_inv(bool if_fast_inv);

    void save_temp_model();
    
    vector<string> SNP_cojo, SNP_PLINK;

    // sumstat: col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
    ArrayXXd sumstat, sumstat_candidate, sumstat_screened, sumstat_backward_new_model;
    

    // X matrix
    vector<bool> X_A1, X_A2;
    vector<double> X_avg, X_std;
    vector<VectorXd> X_candidate, X_backward_new_model;
    
    // necessary information during calculation
    MatrixXd r, R_inv_pre, R_inv_post, R_pre, R_post;
    ArrayXd conditional_beta, beta, beta_var;
    VectorXd X_temp_vec, r_temp_vec;

    // final output
    ArrayXd output_b, output_se2;
    
    double Vp, R2, previous_R2;
    int indi_num;
};


class TCOJO
{
public:
    static long commonSNP_total_num;
    static map<string, int> commonSNP_index_map;
    static vector<int> final_commonSNP_index;
    static vector<string> final_commonSNP;

    static vector<string> A1_ref, A2_ref;
    static vector<int> SNP_pos_ref;

public:
    void read_files_two_cohorts(string cojoFile1, string PLINK1, string cojoFile2, string PLINK2);
    void read_files_one_cohort(string cojoFile, string PLINK);

    void initialize_matrices(Cohort &c);
    void initialize_MDISA(Cohort &c);
    void initialize_backward_selection(Cohort &c, const ArrayXd &pJ);

    void inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
        const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge);
    
    void calc_inner_product_with_candidate(Cohort &c, int single_index);
    void calc_inner_product_with_screened(Cohort &c, int single_index);

    void remove_new_colinear_SNP(bool cohort1_only=false, bool cohort2_only=false);
    void adjust_SNP_according_to_backward_selection(const ArrayXd &pJ, bool cohort1_only=false, bool cohort2_only=false);

    void main_loop(string savename);
    void MDISA(Cohort &c);

    void initialize_hyperparameters(int argc, char** argv);
    void save_results_main_loop(string filepath);
    void save_results_DISA(Cohort &c, string filepath);
    void show_tips_and_exit();

public:
    Cohort c1, c2;

    int max_SNP_index;
    int fixed_candidate_SNP_num=0;

    vector<int> candidate_SNP, screened_SNP, excluded_SNP, backward_removed_SNP;

    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    ArrayXXd sumstat_merge, sumstat_new_model_joint;

    
// hyperparameters for users to predefine and adjust
public: 
    double threshold = 5e-8;
    double colinear_threshold = 0.9;
    double colinear_threshold_sqrt, iter_colinear_threshold;

    double R2_incremental_threshold = 0.0;
    double R2_incremental_threshold_backwards = -0.5;
    
    long window_size = 1500000;
    int max_iter_num = 10000;

    bool if_fast_inv = true;
};

#endif