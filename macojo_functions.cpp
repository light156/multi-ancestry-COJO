#include "macojo.h"


void Cohort::get_vector_from_bed_matrix(int index, ArrayXd &vec)
{
    vec.setZero(indi_num);
    
    const uint64_t words_per_snp = (indi_num + 63) / 64;
    const uint64_t base_word_index = uint64_t(index) * words_per_snp;
    
    const double SNP_avg = X_avg(index);
    int i = 0;

    for (uint64_t w = 0; w < words_per_snp; w++) {
        uint64_t wordA1 = X_A1[base_word_index + w];
        uint64_t wordA2 = X_A2[base_word_index + w];
        int bits_this_word = min(64, indi_num - i);

        for (int b = 0; b < bits_this_word; b++, i++) {
            bool A1 = (wordA1 >> b) & 1ULL;
            bool A2 = (wordA2 >> b) & 1ULL;

            double geno  = double(A1) + double(A2);
            double valid = double((!A1) || A2); // 1.0 if valid, 0.0 otherwise

            if (params.if_keep_NA)
                vec(i) = valid * geno + (1.0 - valid) * -9;
            else
                vec(i) = (geno - SNP_avg) * valid;
        }
    }

    if (!params.if_keep_NA)
        vec /= sqrt(X_square(index));
}


void Cohort::calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int single_index)
{      
    r_temp_vec.setZero(SNP_list.size());
    int temp_index = shared.final_commonSNP_index[single_index];

    if (params.if_LD_mode) {
        for (int j = 0; j < SNP_list.size(); j++) {
            int sweep_index = shared.final_commonSNP_index[SNP_list[j]];
            if (shared.chr_ref[sweep_index] == shared.chr_ref[temp_index] &&
                abs(shared.SNP_pos_ref[sweep_index] - shared.SNP_pos_ref[temp_index]) < params.window_size)
                r_temp_vec(j) = LD_matrix(sweep_index, temp_index);
        }
        return;
    }

    get_vector_from_bed_matrix(temp_index, X_temp_vec);
    
    #pragma omp parallel
    {
        ArrayXd vec_in_list(indi_num);

        #pragma omp for schedule(dynamic)
        for (int j = 0; j < SNP_list.size(); ++j) {
            int sweep_index = shared.final_commonSNP_index[SNP_list[j]];

            if (shared.chr_ref[sweep_index] == shared.chr_ref[temp_index] &&
                abs(shared.SNP_pos_ref[sweep_index] - shared.SNP_pos_ref[temp_index]) < params.window_size) {
                get_vector_from_bed_matrix(sweep_index, vec_in_list);
                r_temp_vec(j) = calc_inner_product(X_temp_vec, vec_in_list, params.if_keep_NA);
            }  
        }
    }
}


void Cohort::save_temp_model() 
{   
    // save current model, irrelevant items may also be saved to avoid complexity
    if (sumstat_candidate.rows() == 1) {
        output_b = sumstat_candidate.col(0);
        output_se2 = sumstat_candidate.col(1);
        R_pre = MatrixXd::Identity(1,1);
        R_inv_pre = MatrixXd::Identity(1,1);
        R_pre_gcta = MatrixXd::Identity(1,1) * sumstat_candidate(0,6); // only used when if_gcta_COJO == true
        R_inv_pre_gcta = MatrixXd::Identity(1,1) / sumstat_candidate(0,6); // only used when if_gcta_COJO == true
        previous_R2 = 0; // only used for if_gcta_COJO == false;
    } else {
        output_b = beta;
        output_se2 = beta_var;
        R_inv_pre = R_inv_post;
        R_pre = R_post; // only used for if_fast_inv == false
        R_inv_pre_gcta = R_inv_post_gcta; // only used when if_gcta_COJO == true
        R_pre_gcta = R_post_gcta; // only used when if_gcta_COJO == true
        previous_R2 = R2; // only used for if_gcta_COJO == false
    }
}


void Cohort::calc_conditional_effects() 
{   
    if (params.if_gcta_COJO) {
        MatrixXd temp1 = sumstat_candidate.col(0) * sumstat_candidate.col(6);
        MatrixXd temp2 = R_inv_pre_gcta * temp1;
        ArrayXd temp3 = r_gcta * temp2;
        conditional_beta = sumstat_screened.col(0) - temp3 / sumstat_screened.col(6);
    } else {
        MatrixXd temp1 = sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5));
        MatrixXd temp2 = R_inv_pre * temp1;
        ArrayXd temp3 = r * temp2;
        conditional_beta = sumstat_screened.col(0) - temp3 / sqrt(sumstat_screened.col(5));
    }
}


// should be called after sumstat_new_model is updated
bool Cohort::calc_joint_effects()
{   
    if (params.if_gcta_COJO) {
        VectorXd temp1 = sumstat_candidate.col(6) * sumstat_candidate.col(0);
        beta = R_inv_post_gcta * temp1;
        beta_var = R_inv_post_gcta.diagonal().array() * Vp;
    } else {
        VectorXd temp1 = sqrt(sumstat_candidate.col(5)) * sumstat_candidate.col(0);
        ArrayXd temp2 = R_inv_post * temp1;
        beta = temp2 / sqrt(sumstat_candidate.col(5));

        double sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
        beta_var = sigma_J_squared * R_inv_post.diagonal().array() / sumstat_candidate.col(6);

        double Neff = median(sumstat_candidate.col(4));
        int M = sumstat_candidate.rows();
        R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    }

    return beta_var.minCoeff() > 1e-30;
}


bool Cohort::calc_R_inv() 
{   
    // if (params.if_fast_inv) 
    //    return calc_R_inverse_fast(R_inv_pre, r_temp_vec, 1.0, params.iter_colinear_threshold, R_inv_post);
    
    ArrayXd scaling_vector;
    if (params.if_keep_NA || params.if_LD_mode || !params.if_gcta_COJO) 
        scaling_vector = ArrayXd::Ones(sumstat_candidate.rows());
    else 
        scaling_vector = sumstat_candidate.col(8);

    return calc_R_inverse_forward(R_pre, r_temp_vec, 1.0, params.iter_colinear_threshold, R_post, R_inv_post, true, scaling_vector);
}


void Cohort::calc_R_inv_gcta(int screened_index) 
{   
    r_temp_vec_gcta = r_gcta.row(screened_index).transpose();
    double lower_right_corner = sumstat_screened(screened_index, 6);

    calc_R_inverse_forward(R_pre_gcta, r_temp_vec_gcta, lower_right_corner, params.iter_colinear_threshold, 
                    R_post_gcta, R_inv_post_gcta, false);

    // if (params.if_fast_inv) 
    //    calc_R_inverse_fast(R_inv_pre_gcta, r_temp_vec_gcta, lower_right_corner, params.iter_colinear_threshold, R_inv_post_gcta);
}


bool Cohort::calc_R_inv_from_SNP_list(const vector<int> &SNP_list, const ArrayXXd &sumstat) 
{   
    if (SNP_list.size() != sumstat.rows()) 
        LOGGER.e(0, "Error: SNP list size not equal to sumstat rows");

    int total_num = SNP_list.size();
    R_post = MatrixXd::Identity(total_num, total_num);

    if (params.if_LD_mode) {
        for (int i = 0; i < total_num; i++) {
            int index_i = shared.final_commonSNP_index[SNP_list[i]];
            for (int j = 0; j < i; j++) {
                int index_j = shared.final_commonSNP_index[SNP_list[j]];
                if (shared.chr_ref[index_i] == shared.chr_ref[index_j] &&
                    abs(shared.SNP_pos_ref[index_i] - shared.SNP_pos_ref[index_j]) < params.window_size) {
                    R_post(i, j) = LD_matrix(index_i, index_j);
                    R_post(j, i) = R_post(i, j);
                }
            }
        }
    } else {
        vector<ArrayXd> X_new_model;
        
        for (int index : SNP_list) {
            get_vector_from_bed_matrix(shared.final_commonSNP_index[index], X_temp_vec);
            X_new_model.push_back(X_temp_vec);
        }

        for (int i = 0; i < total_num; i++) {
            int index_i = shared.final_commonSNP_index[SNP_list[i]];
            for (int j = 0; j < i; j++) {
                int index_j = shared.final_commonSNP_index[SNP_list[j]];
                if (shared.chr_ref[index_i] == shared.chr_ref[index_j] &&
                    abs(shared.SNP_pos_ref[index_i] - shared.SNP_pos_ref[index_j]) < params.window_size) {
                    R_post(i, j) = calc_inner_product(X_new_model[i], X_new_model[j], params.if_keep_NA);
                    R_post(j, i) = R_post(i, j);
                }
            }
        }
    } 
    
    R_inv_post = R_post.ldlt().solve(MatrixXd::Identity(total_num, total_num));

    if (params.if_gcta_COJO) {
        R_post_gcta.setZero(total_num, total_num);

        for (int i = 0; i < total_num; i++) {
            R_post_gcta(i, i) = sumstat(i, 6);

            for (int j = 0; j < i; j++) {
                R_post_gcta(i, j) = R_post(i, j) * min(sumstat(i, 4), sumstat(j, 4)) * \
                    sqrt(sumstat(i, 5) * sumstat(j, 5));
                R_post_gcta(j, i) = R_post_gcta(i, j);
            }
        }

        R_inv_post_gcta = R_post_gcta.ldlt().solve(MatrixXd::Identity(total_num, total_num));
    }

    return R_inv_post.cwiseAbs().maxCoeff() < params.iter_colinear_threshold;
}


void MACOJO::accept_SNP_as_candidate(int screened_index) 
{   
    int candidate_index = screened_SNP[screened_index];
    candidate_SNP.push_back(candidate_index);
    screened_SNP.erase(screened_SNP.begin() + screened_index);

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];
        remove_row(c.sumstat_screened, screened_index);
        remove_row(c.r, screened_index);
        if (params.if_gcta_COJO) 
            remove_row(c.r_gcta, screened_index);
        
        if (screened_SNP.size() == 0) continue;

        c.calc_inner_product_with_SNP_list(screened_SNP, candidate_index);
        append_column(c.r, c.r_temp_vec);

        if (params.if_gcta_COJO) {
            double N_new = c.sumstat_candidate(c.sumstat_candidate.rows()-1, 4);
            double V_new = c.sumstat_candidate(c.sumstat_candidate.rows()-1, 5);
            c.r_temp_vec_gcta = c.r_temp_vec.array() * c.sumstat_screened.col(4).min(N_new) * \
                sqrt(c.sumstat_screened.col(5)) * sqrt(V_new);
            append_column(c.r_gcta, c.r_temp_vec_gcta);
        }
    }
}


void MACOJO::remove_SNP_from_screened(int screened_index) 
{   
    int removed_index = screened_SNP[screened_index];
    removed_SNP.push_back(removed_index);
    screened_SNP.erase(screened_SNP.begin() + screened_index);

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];
        append_row(c.sumstat_removed, c.sumstat_screened.row(screened_index));
        remove_row(c.sumstat_screened, screened_index);
        remove_row(c.r, screened_index);

        if (params.if_gcta_COJO)
            remove_row(c.r_gcta, screened_index);
    }
}


void MACOJO::remove_SNP_from_candidate(int candidate_index) 
{   
    int removed_index = candidate_SNP[candidate_index];
    removed_SNP.push_back(removed_index);
    candidate_SNP.erase(candidate_SNP.begin() + candidate_index);

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];
        append_row(c.sumstat_removed, c.sumstat_candidate.row(candidate_index));
        remove_row(c.sumstat_candidate, candidate_index);
        remove_column(c.r, candidate_index);

        if (params.if_gcta_COJO)
            remove_column(c.r_gcta, candidate_index);
    }
}


void MACOJO::initialize_MDISA() 
{   
    for (int n : current_calculation_list) {
        auto &c = cohorts[n];
        
        for (int i = 0; i < removed_SNP.size(); i++) {
            append_row(c.sumstat_screened, c.sumstat_removed.row(i));

            c.calc_inner_product_with_SNP_list(candidate_SNP, removed_SNP[i]);
            append_row(c.r, c.r_temp_vec.transpose());

            if (params.if_gcta_COJO) {
                double N_new = c.sumstat_removed(i, 4);
                double V_new = c.sumstat_removed(i, 5);
                c.r_temp_vec_gcta = c.r_temp_vec.array() * c.sumstat_candidate.col(4).min(N_new) * \
                    sqrt(c.sumstat_candidate.col(5)) * sqrt(V_new);

                append_row(c.r_gcta, c.r_temp_vec_gcta.transpose());
            }
        }

        c.sumstat_removed.resize(0, 0);
    }

    screened_SNP.insert(screened_SNP.end(), removed_SNP.begin(), removed_SNP.end());
    removed_SNP.clear();
}


void MACOJO::inverse_var_meta_init() 
{
    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        abs_zC = abs(c.sumstat_screened.col(0)) / sqrt(c.sumstat_screened.col(1));
        return;
    }

    ArrayXd num = ArrayXd::Zero(screened_SNP.size()), denom = ArrayXd::Zero(screened_SNP.size());

    for (int n : current_calculation_list) {
        auto &c = cohorts[current_calculation_list[n]];
        num += c.sumstat_screened.col(0) / c.sumstat_screened.col(1);
        denom += 1 / c.sumstat_screened.col(1);
    }
    
    abs_zC = abs(num) / sqrt(denom);
}


void MACOJO::inverse_var_meta_conditional() 
{
    ArrayXd num = ArrayXd::Zero(screened_SNP.size()), denom = ArrayXd::Zero(screened_SNP.size());

    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        num = c.conditional_beta;
        if (params.if_gcta_COJO)
            denom = c.Vp / c.sumstat_screened.col(6);
        else
            denom = c.sumstat_screened.col(1);
    } else {
        ArrayXd temp_se2;
        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]]; 
            if (params.if_gcta_COJO) 
                temp_se2 = c.Vp / c.sumstat_screened.col(6);
            else
                temp_se2 = c.sumstat_screened.col(1);

            num += c.conditional_beta / temp_se2;
            denom += 1.0 / temp_se2;
        }   
    }

    abs_zC = abs(num) / sqrt(denom);
}


void MACOJO::inverse_var_meta_joint() 
{
    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        bJ = c.beta;
        se2J = c.beta_var;
    } else {
        bJ.setZero(cohorts[current_calculation_list[0]].beta.rows());
        se2J.setZero(cohorts[current_calculation_list[0]].beta.rows());

        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]];
            bJ += c.beta / c.beta_var;
            se2J += 1 / c.beta_var;
        }

        se2J = 1 / se2J;
        bJ = bJ * se2J;
    }   

    abs_zJ = abs(bJ) / sqrt(se2J);
}


void MACOJO::output_results_to_file(string filepath) 
{   
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    jmaCOJO << "Chr"
            << "\t" << "SNP" 
            << "\t" << "bp" 
            << "\t" << "A1" 
            << "\t" << "A2";

    if (current_calculation_list.size() == 1) {
        jmaCOJO << "\t" << "freq"
                << "\t" << "b"
                << "\t" << "se"
                << "\t" << "p"
                << "\t" << "n"
                << "\t" << "bJ"
                << "\t" << "bJ_se"
                << "\t" << "pJ";
    } else {
        for (int n : current_calculation_list) {
            jmaCOJO << "\t" << "freq." << to_string(n+1) 
                    << "\t" << "b." << to_string(n+1) 
                    << "\t" << "se." << to_string(n+1) 
                    << "\t" << "p." << to_string(n+1) 
                    << "\t" << "n." << to_string(n+1) 
                    << "\t" << "bJ." << to_string(n+1) 
                    << "\t" << "bJ_se." << to_string(n+1) 
                    << "\t" << "pJ." << to_string(n+1);
        }
    
        jmaCOJO << "\t" << "bJ.ma"  
                << "\t" << "bJ_se.ma" 
                << "\t" << "pJ.ma";
                
        inverse_var_meta_joint(); 
    }
    
    jmaCOJO << "\n";

    int index = 0;
    map<int, int> SNP_ref_order_pair;
 
    // Inserting element in pair vector to keep track of previous indexes
    for (auto SNP_index : candidate_SNP) {
        SNP_ref_order_pair.insert(make_pair(shared.SNP_pos_ref[shared.final_commonSNP_index[SNP_index]], index));
        index++;
    }
    
    jmaCOJO.precision(12);

    for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) {
        index = iter->second;
        int ref_index = shared.final_commonSNP_index[candidate_SNP[index]];
        jmaCOJO << shared.chr_ref[ref_index]
                << "\t" << shared.final_commonSNP[candidate_SNP[index]] 
                << "\t" << iter->first 
                << "\t" << shared.A1_ref[ref_index] 
                << "\t" << shared.A2_ref[ref_index];

        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            
            // cols 0:b, 1:se2, 2:p, 3:freq, 4:N, 5: bJ, 6:seJ, 7:pJ
            jmaCOJO << "\t" << c.sumstat_candidate(index, 3) 
                    << "\t" << c.sumstat_candidate(index, 0) 
                    << "\t" << sqrt(c.sumstat_candidate(index, 1)) 
                    << "\t" << c.sumstat_candidate(index, 2) 
                    << "\t" << c.sumstat_candidate(index, 4) 

                    << "\t" << c.output_b(index) 
                    << "\t" << sqrt(c.output_se2(index)) 
                    << "\t" << erfc(abs(c.output_b(index)) / sqrt(c.output_se2(index)) / sqrt(2));
        }

        if (current_calculation_list.size() > 1) {
            jmaCOJO << "\t" << bJ(index) 
                    << "\t" << se2J(index)
                    << "\t" << erfc(abs_zJ(index) / sqrt(2));
        } 
        
        jmaCOJO << "\n";
    }

    jmaCOJO.clear();
    jmaCOJO.close();
    LOGGER.i(0, "Results saved into [" + filepath + "]\n");
}


void MACOJO::output_user_hyperparameters() 
{   
    LOGGER << endl << "=========== MACOJO CONFIGURATION ===========" << endl;
    LOGGER << "p-value Threshold: " << params.threshold << endl
            << "Colinearity threshold: " << params.colinear_threshold << endl
            << "R2 incremental threshold: " << params.R2_incremental_threshold << endl
            << "R2 incremental threshold backwards: " << params.R2_incremental_threshold_backwards << endl
            << "SNP position window (+/-): " << params.window_mb << "Mb" << endl
            << "SNP frequency threshold: " << params.freq_threshold << endl
            << (cohorts.size() > 1 ? 
                (params.if_freq_mode_and ? "SNP frequency mode: AND\n" : "SNP frequency mode: OR\n") : "")
            << (params.if_LD_mode ? "Use .ld files for calculation\n" : "")
            << (params.if_MDISA ? "Run MDISA after MACOJO\n" : "")
            << (params.if_keep_NA ? "Do not fill NA with mean genotype values\n" : "")
            << (params.if_gcta_COJO ? "Use GCTA-COJO model selection criteria\n" : "")
            << "===========================================" << endl << endl;
}
