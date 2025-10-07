#pragma once
#include "BitArray.h"
#include "LD.h"
#include "Logger.h"
#include <omp.h>
#include <map>
#include <vector>
#include <list>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/SpecialFunctions>
#include <iomanip>

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
    void read_sumstat(string cojofile);
    void read_PLINK(string PLINKfile, bool is_ref_cohort);
    void read_PLINK_LD(string PLINKfile, bool is_ref_cohort);
    void skim_fam(string famFile);

    void get_vector_from_bed_matrix(int index, VectorXd &vec);
    void calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int single_index);

    void calc_conditional_effects();
    bool calc_joint_effects(const ArrayXXd &sumstat_temp);    
    void calc_Vp();
    void calc_R_inv();
    void calc_R_inv_from_SNP_list(const vector<int> &SNP_list);

public:
    // sumstat: col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
    ArrayXXd sumstat, sumstat_candidate, sumstat_screened, sumstat_backward_new_model;
    
    // X or LD matrix, depend on which mode is used
    BitArray X_A1, X_A2;
    vector<double> X_avg, X_std;
    LDPacked LD_matrix;

    // necessary information during calculation
    MatrixXd r, R_inv_post, R_post, R_inv_pre = MatrixXd::Identity(1,1), R_pre = MatrixXd::Identity(1,1);
    ArrayXd conditional_beta, beta, beta_var;
    VectorXd X_temp_vec, r_temp_vec;

    // final output
    ArrayXd output_b, output_se2;
    
    double Vp, R2, previous_R2 = 0.0;
    uint64_t indi_num;

public:
    int window_size;
    bool if_fast_inv;
    bool if_LD_mode;
    bool if_fill_NA;
    double iter_colinear_threshold;
};


class MACOJO
{
public:
    static uint64_t commonSNP_total_num;
    static map<string, int> commonSNP_index_map;
    static vector<string> A1_ref, A2_ref;
    static vector<int> SNP_pos_ref;
    static vector<int> final_commonSNP_index;
    vector<string> final_commonSNP;

public:
    void read_user_hyperparameters(int argc, char** argv);
    void read_cojo_PLINK_files(char** filenames, int cohort_num);
    void read_SNP_only(string filename, vector<string> &SNP_list, bool if_sumstat=false, bool if_bim=false);
    void set_reference_from_bim(string PLINKfile);
    void entry_function(string savename);

    void initialize_main_loop();
    void main_loop();
    void inverse_var_meta(bool if_conditional);
    void inverse_var_meta_joint();
    void accept_SNP_as_candidate(int screened_index);
    void adjust_SNP_according_to_backward_selection();
    void save_temp_model();
    
    void initialize_MDISA_from_MACOJO(Cohort &c);
    void output_results_to_file(string filepath);

    void show_tips();

public:
    vector<Cohort> cohorts;
    vector<int> current_calculation_list;
    vector<int> candidate_SNP, screened_SNP, backward_removed_SNP;
    
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    ArrayXXd sumstat_merge, sumstat_new_model_joint;
    
// hyperparameters for users to predefine and adjust
public: 
    double threshold = 5e-8;
    double colinear_threshold = 0.9;
    double freq_threshold = 0.01;
    bool if_freq_mode_or = true;
    double R2_incremental_threshold = -1;
    double R2_incremental_threshold_backwards = -1;
    int window_size = 1e7;
    int max_iter_num = 10000;

    bool if_fast_inv = true;
    bool if_MDISA = true;
    bool if_cojo_joint = false;
    bool if_LD_mode = false;
    bool if_fill_NA = false;

public:
    vector<string> all_SNP, fixed_candidate_SNP;
    int fixed_candidate_SNP_num=0;
};
