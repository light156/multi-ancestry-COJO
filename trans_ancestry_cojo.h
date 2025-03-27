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


class Cohort {
    
public:
    void read_sumstat(string cojofile);
    void read_PLINK(string PLINKfile, bool is_ref_cohort);
    void generate_sumstat_and_X();

    void calc_inner_product(const vector<int> &index_list, int single_index, int window_size);
    void calc_conditional_effects();
    bool calc_joint_effects(const ArrayXXd &sumstat, bool flag, double iter_colinear_threshold=0);
    void calc_R_inv_fast();
    
    static map<string, int> sumstat_commonSNP_index;
    static vector<string> commonSNP;
    static vector<string> A1_ref, A2_ref;
    static vector<int> SNP_pos_ref;

    MatrixXd X;
    ArrayXXd sumstat, sumstat_candidate, sumstat_screened, sumstat_backward;
    
    MatrixXd r, R_inv_pre, R_inv_post;
    VectorXd r_temp_vec;
    ArrayXd conditional_beta, beta, beta_var;

    // final output
    ArrayXd output_b, output_se2;
    
    double Vp, R2;
    int indi_num;

public:
    map<string, int> SNP_index;
    vector<string> A1, A2;
    vector<double> freq, b, se, p, N;
    vector<vector<double>> genotype;
    vector<string> SNP, included_SNP_PLINK;
};


class TransAncestryCOJO {

public:
    void read_files(string cojoFile1, string cojoFile2, string PLINK1, string PLINK2);

    void main_loop(string savename);
    void save_results_main_loop(string filepath);

    void MDISA(Cohort &c);
    void save_results_DISA(Cohort &c, string filepath);

// Step 2: main loop
public:
    void initialize_hyperparameters();

    void inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
        const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge);

    void initialize_main_loop(Cohort &c); 
    void initialize_MDISA(Cohort &c);
    void initialize_backward_selection(Cohort &c);
    
    void remove_new_colinear_SNP(bool cohort1_only=false, bool cohort2_only=false);
    void adjust_SNP_according_to_backward_selection(bool cohort1_only=false, bool cohort2_only=false);

public:
    Cohort c1, c2;

    int max_SNP_index;
    int fixed_candidate_SNP_num=0;
    vector<int> candidate_SNP, screened_SNP, excluded_SNP, backward_removed_SNP;

    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    ArrayXXd sumstat_merge, sumstat_new_model_joint;

    // final output 
    ArrayXd output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2;
    

// hyperparameters for users to predefine and adjust
public: 
    double threshold;
    double colinear_threshold, colinear_threshold_sqrt, iter_colinear_threshold;

    double R2_incremental_threshold;
    double R2_incremental_threshold_backwards;
    
    int window_size;
    int max_iter_num;
};

#endif
