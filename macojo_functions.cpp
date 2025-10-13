#include "macojo.h"


void Cohort::get_vector_from_bed_matrix(int index, ArrayXd &vec)
{
    vec.setZero(indi_num);
    uint64_t X_start_index = uint64_t(index) * indi_num;
    uint64_t start_word = X_start_index >> 6;
    uint64_t start_bit  = X_start_index & 63;

    double SNP_avg = X_avg(index);
    int i = 0;

    while (i < indi_num) {
        uint64_t wordA1 = X_A1.data_[start_word];
        uint64_t wordA2 = X_A2.data_[start_word];

        for (uint64_t b = start_bit; b < 64 && i < indi_num; b++, i++) {
            bool A1 = (wordA1 >> b) & 1ULL;
            bool A2 = (wordA2 >> b) & 1ULL;

            if (if_keep_NA)
                vec(i) = (!A1 || A2) ? double(A1) + double(A2) : -9;
            else
                vec(i) = (!A1 || A2) ? double(A1) + double(A2) - SNP_avg : 0;
        }

        start_word++;
        start_bit = 0;
    }

    if (!if_keep_NA)
        vec /= sqrt(X_square(index));
}


void Cohort::calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int single_index)
{      
    r_temp_vec.setZero(SNP_list.size());
    int temp_index = MACOJO::final_commonSNP_index[single_index];
    int temp_SNP_pos = MACOJO::SNP_pos_ref[temp_index];

    if (if_LD_mode) {
        for (int j = 0; j < SNP_list.size(); j++) {
            if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]] - temp_SNP_pos) < window_size)
                r_temp_vec(j) = LD_matrix(MACOJO::final_commonSNP_index[SNP_list[j]], temp_index);
        }
    } else {
        get_vector_from_bed_matrix(temp_index, X_temp_vec);

        #pragma omp parallel for
        for (int j = 0; j < SNP_list.size(); j++) {
            if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]] - temp_SNP_pos) < window_size) {
                ArrayXd vec_in_list;
                get_vector_from_bed_matrix(MACOJO::final_commonSNP_index[SNP_list[j]], vec_in_list);
                r_temp_vec(j) = calc_inner_product(X_temp_vec, vec_in_list, if_keep_NA);
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
        previous_R2 = R2;
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
    if (if_gcta_COJO) {
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
    if (if_gcta_COJO) {
        VectorXd temp1 = sumstat_new_model.col(6) * sumstat_new_model.col(0);
        beta = R_inv_post_gcta * temp1;
        beta_var = R_inv_post_gcta.diagonal().array() * Vp;
    } else {
        VectorXd temp1 = sqrt(sumstat_new_model.col(5)) * sumstat_new_model.col(0);
        ArrayXd temp2 = R_inv_post * temp1;
        beta = temp2 / sqrt(sumstat_new_model.col(5));

        double sigma_J_squared = Vp - (sumstat_new_model.col(0) * sumstat_new_model.col(5) * beta).sum();
        beta_var = sigma_J_squared * R_inv_post.diagonal().array() / sumstat_new_model.col(6);

        double Neff = median(sumstat_new_model.col(4));
        int M = sumstat_new_model.rows();
        R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    }

    return beta_var.minCoeff() > 1e-30;
}


bool Cohort::calc_R_inv() 
{   
    if (if_fast_inv) 
        return calc_block_inverse_fast(R_inv_pre, r_temp_vec, 1.0, iter_colinear_threshold, R_inv_post);
    else {
        ArrayXd scaling_vector = sumstat_candidate.col(8);
        return calc_block_inverse_exact(R_pre, r_temp_vec, 1.0, iter_colinear_threshold, R_post, R_inv_post, true, &scaling_vector);
    }
}


void Cohort::calc_R_inv_gcta(int screened_index) 
{   
    r_temp_vec_gcta = r_gcta.row(screened_index).transpose();
    double lower_right_corner = sumstat_screened(screened_index, 6);

    if (if_fast_inv) 
        calc_block_inverse_fast(R_inv_pre_gcta, r_temp_vec_gcta, lower_right_corner, iter_colinear_threshold, R_inv_post_gcta);
    else
        calc_block_inverse_exact(R_pre_gcta, r_temp_vec_gcta, lower_right_corner, iter_colinear_threshold, R_post_gcta, R_inv_post_gcta, false);
}


void Cohort::calc_R_inv_backward(int remove_index) 
{   
    if (if_fast_inv) {
        calc_block_inverse_fast(R_inv_pre, remove_index, R_inv_post);
        if (if_gcta_COJO)
            calc_block_inverse_fast(R_inv_pre_gcta, remove_index, R_inv_post_gcta);
    } else {
        calc_block_inverse_exact(R_pre, remove_index, R_post, R_inv_post);
        if (if_gcta_COJO)
            calc_block_inverse_exact(R_pre_gcta, remove_index, R_post_gcta, R_inv_post_gcta);
    }
}


bool Cohort::calc_R_inv_from_SNP_list(const vector<int> &SNP_list, const ArrayXXd &sumstat) 
{   
    if (SNP_list.size() != sumstat_new_model.rows()) 
        LOGGER.e(0, "Error: SNP list size not equal to sumstat rows");

    int total_num = SNP_list.size();
    R_post = MatrixXd::Identity(total_num, total_num);

    if (if_LD_mode) {
        for (int i = 0; i < total_num; i++) {
            for (int j = 0; j < i; j++) {
                if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[i]]] - 
                        MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]]) < window_size) {
                    R_post(i, j) = LD_matrix(MACOJO::final_commonSNP_index[SNP_list[i]], MACOJO::final_commonSNP_index[SNP_list[j]]);
                    R_post(j, i) = R_post(i, j);
                }
            }
        }
    } else {
        vector<ArrayXd> X_new_model;
        
        for (int index : SNP_list) {
            get_vector_from_bed_matrix(MACOJO::final_commonSNP_index[index], X_temp_vec);
            X_new_model.push_back(X_temp_vec);
        }

        for (int i = 0; i < total_num; i++) {
            for (int j = 0; j < i; j++) {
                if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[i]]] - 
                        MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]]) < window_size) {             
                    R_post(i, j) = calc_inner_product(X_new_model[i], X_new_model[j], if_keep_NA);
                    R_post(j, i) = R_post(i, j);
                }
            }
        }
    } 
    
    R_inv_post = R_post.ldlt().solve(MatrixXd::Identity(total_num, total_num));

    if (if_gcta_COJO) {
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

    return R_inv_post.cwiseAbs().maxCoeff() < iter_colinear_threshold;
}


void MACOJO::accept_SNP_as_candidate(int screened_index) 
{   
    int candidate_index = screened_SNP[screened_index];
    candidate_SNP.push_back(candidate_index);
    screened_SNP.erase(screened_SNP.begin() + screened_index);

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];
        append_row(c.sumstat_candidate, c.sumstat_screened.row(screened_index));
        remove_row(c.sumstat_screened, screened_index);

        c.calc_inner_product_with_SNP_list(screened_SNP, candidate_index);
        remove_row(c.r, screened_index);
        append_column(c.r, c.r_temp_vec);

        if (if_gcta_COJO) {
            double N_new = c.sumstat_candidate(c.sumstat_candidate.rows()-1, 4);
            double V_new = c.sumstat_candidate(c.sumstat_candidate.rows()-1, 5);
            c.r_temp_vec_gcta = c.r_temp_vec.array() * c.sumstat_screened.col(4).min(N_new) * \
                sqrt(c.sumstat_screened.col(5)) * sqrt(V_new);

            remove_row(c.r_gcta, screened_index);
            append_column(c.r_gcta, c.r_temp_vec_gcta);
        }
    }
}


void MACOJO::initialize_MDISA() 
{   
    for (int i = removed_SNP.size()-1; i >= 0; i--) {
        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            append_row(c.sumstat_screened, c.sumstat_removed.row(i));
            remove_row(c.sumstat_removed, i);

            c.calc_inner_product_with_SNP_list(candidate_SNP, removed_SNP[i]);
            append_row(c.r, c.r_temp_vec.transpose());

            if (if_gcta_COJO) {
                double N_new = c.sumstat_removed(i, 4);
                double V_new = c.sumstat_removed(i, 5);
                c.r_temp_vec_gcta = c.r_temp_vec.array() * c.sumstat_candidate.col(4).min(N_new) * \
                    sqrt(c.sumstat_candidate.col(5)) * sqrt(V_new);

                append_row(c.r_gcta, c.r_temp_vec_gcta.transpose());
            }
        }

        screened_SNP.push_back(removed_SNP[i]);
        removed_SNP.erase(removed_SNP.begin() + i);
    }
}


void MACOJO::inverse_var_meta_init() 
{
    // merge: 2:Zabs, 3:p
    sumstat_merge.setZero(screened_SNP.size(), 4);

    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        sumstat_merge.col(0) = c.sumstat_screened.col(0);
        sumstat_merge.col(1) = c.sumstat_screened.col(1);
    } else {
        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]];
            sumstat_merge.col(0) += c.sumstat_screened.col(0) / c.sumstat_screened.col(1);
            sumstat_merge.col(1) += 1 / c.sumstat_screened.col(1);
        }
    }

    sumstat_merge.col(2) = abs(sumstat_merge.col(0)) / sqrt(sumstat_merge.col(1));
    sumstat_merge.col(3) = erfc(sumstat_merge.col(2) / sqrt(2));
}


void MACOJO::inverse_var_meta_conditional() 
{
    // merge: 2:Zabs, 3:p
    sumstat_merge.setZero(screened_SNP.size(), 4);

    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        sumstat_merge.col(0) = c.conditional_beta;
        if (if_gcta_COJO)
            sumstat_merge.col(1) = c.Vp / c.sumstat_screened.col(6);
        else
            sumstat_merge.col(1) = c.sumstat_screened.col(1);

    } else {
        ArrayXd temp_se2;

        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]]; 
            if (if_gcta_COJO) 
                temp_se2 = c.Vp / c.sumstat_screened.col(6);
            else
                temp_se2 = c.sumstat_screened.col(1);

            sumstat_merge.col(0) += c.conditional_beta / temp_se2;
            sumstat_merge.col(1) += 1.0 / temp_se2;
        }   
    }

    sumstat_merge.col(2) = abs(sumstat_merge.col(0)) / sqrt(sumstat_merge.col(1));
    sumstat_merge.col(3) = erfc(sumstat_merge.col(2) / sqrt(2));
}


void MACOJO::inverse_var_meta_joint() 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    sumstat_merge_new_model.setZero(cohorts[current_calculation_list[0]].beta.rows(), 4);
    
    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        sumstat_merge_new_model.col(0) = c.beta;
        sumstat_merge_new_model.col(1) = c.beta_var;
    } else {
        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]];
            sumstat_merge_new_model.col(0) += c.beta / c.beta_var;
            sumstat_merge_new_model.col(1) += 1 / c.beta_var;
        }

        sumstat_merge_new_model.col(1) = 1 / sumstat_merge_new_model.col(1);
        sumstat_merge_new_model.col(0) = sumstat_merge_new_model.col(0) * sumstat_merge_new_model.col(1);
    }

    sumstat_merge_new_model.col(2) = abs(sumstat_merge_new_model.col(0)) / sqrt(sumstat_merge_new_model.col(1));
    sumstat_merge_new_model.col(3) = erfc(sumstat_merge_new_model.col(2) / sqrt(2));
}


void MACOJO::output_results_to_file(string filepath) 
{   
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    jmaCOJO << "SNP" 
            << "\t" << "bp" 
            << "\t" << "A1" 
            << "\t" << "A2";

    for (int n : current_calculation_list)
        jmaCOJO << "\t" << "freq." << to_string(n+1) 
                << "\t" << "b." << to_string(n+1) 
                << "\t" << "se." << to_string(n+1) 
                << "\t" << "p." << to_string(n+1) 
                << "\t" << "N." << to_string(n+1) 
                << "\t" << "bJ." << to_string(n+1) 
                << "\t" << "seJ." << to_string(n+1) 
                << "\t" << "pJ." << to_string(n+1);

    ArrayXd bJma, seJma, pJma;

    if (current_calculation_list.size() > 1) {
        jmaCOJO << "\t" << "bJ.ma"  
                << "\t" << "seJ.ma" 
                << "\t" << "pJ.ma";
                
        // calculate output sumstats    
        sumstat_merge_new_model.setZero(candidate_SNP.size(), 4);

        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            sumstat_merge_new_model.col(0) += c.output_b / c.output_se2;
            sumstat_merge_new_model.col(1) += 1 / c.output_se2;
        }
        
        bJma = sumstat_merge_new_model.col(0) / sqrt(sumstat_merge_new_model.col(1));
        seJma = 1 / sqrt(sumstat_merge_new_model.col(1));
        pJma = erfc(abs(bJma) / seJma / sqrt(2));
    } 
    
    jmaCOJO << "\n";

    int index = 0;
    map<int, int> SNP_ref_order_pair;
 
    // Inserting element in pair vector to keep track of previous indexes
    for (auto SNP_index : candidate_SNP) {
        SNP_ref_order_pair.insert(make_pair(SNP_pos_ref[final_commonSNP_index[SNP_index]], index));
        index++;
    }
    
    jmaCOJO.precision(12);

    for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) {
        index = iter->second;
        jmaCOJO << final_commonSNP[index] 
                << "\t" << iter->first 
                << "\t" << A1_ref[final_commonSNP_index[candidate_SNP[index]]] 
                << "\t" << A2_ref[final_commonSNP_index[candidate_SNP[index]]];

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
            jmaCOJO << "\t" << bJma(index) 
                    << "\t" << seJma(index)
                    << "\t" << pJma(index);
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
    LOGGER << "Threshold: 5e-8" << endl
            << "Colinearity threshold: " << colinear_threshold << endl
            << "R2 incremental threshold: " << R2_incremental_threshold << endl
            << "R2 incremental threshold backwards: " << R2_incremental_threshold_backwards << endl
            << "SNP position window (+/-): " << window_mb << "Mb" << endl
            << "SNP frequency threshold: " << freq_threshold << endl
            << "SNP frequency mode: " << (if_freq_mode_and ? "AND" : "OR") << endl
            << (if_LD_mode ? "Use .ld files for calculation\n" : "")
            << (if_skip_MDISA ? "Do not run MDISA after MACOJO\n" : "")
            << (if_keep_NA ? "Do not fill NA with mean genotype values\n" : "")
            << (if_gcta_COJO ? "Use GCTA-COJO model selection criteria\n" : "")
            << "===========================================" << endl << endl;

    for (auto &c : cohorts) {
        c.window_size = window_mb > 0 ? window_mb * 1e6 : INT_MAX;
        c.if_fast_inv = if_fast_inv;
        c.if_LD_mode = if_LD_mode;
        c.if_keep_NA = if_keep_NA;
        c.if_gcta_COJO = if_gcta_COJO;
        c.iter_colinear_threshold = 1 / (1 - colinear_threshold);
    }
}
