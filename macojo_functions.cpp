#include "macojo.h"


void Cohort::get_vector_from_bed_matrix(int index, VectorXd &vec)
{
    vec.setZero(indi_num);
    uint64_t X_start_index = uint64_t(index) * indi_num;

    double mean = X_avg[index], stdv = X_std[index];
    /*
    for (int i = 0; i < indi_num; i++) {
        bool A1 = X_A1.get(X_start_index + i), A2 = X_A2.get(X_start_index + i);
        if (!A1 || A2) 
            vec(i) = (double(A1) + double(A2) - mean) / stdv;
    }
    */

    uint64_t start_word = X_start_index >> 6;
    uint64_t start_bit  = X_start_index & 63;

    int i = 0;

    while (i < indi_num) {
        uint64_t wordA1 = X_A1.data_[start_word];
        uint64_t wordA2 = X_A2.data_[start_word];

        for (uint64_t b = start_bit; b < 64 && i < indi_num; b++, i++) {
            bool A1 = (wordA1 >> b) & 1ULL;
            bool A2 = (wordA2 >> b) & 1ULL;

            if (!A1 || A2)
                vec(i) = (double(A1) + double(A2) - mean) / stdv;
        }

        start_word++;
        start_bit = 0;
    }
}


void Cohort::calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int single_index, int window_size, bool if_LD_mode)
{   
    r_temp_vec.setZero(SNP_list.size());
    int temp_index = MACOJO::final_commonSNP_index[single_index];
    int temp_SNP_pos = MACOJO::SNP_pos_ref[temp_index];

    if (if_LD_mode) {
        for (int j = 0; j < SNP_list.size(); j++) {
            if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]] - temp_SNP_pos) <= window_size)
                r_temp_vec(j) = LD_matrix(MACOJO::final_commonSNP_index[SNP_list[j]], temp_index);
        }
    } else {
        get_vector_from_bed_matrix(temp_index, X_temp_vec);

        #pragma omp parallel for
        for (int j = 0; j < SNP_list.size(); j++) {
            if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]] - temp_SNP_pos) <= window_size) {
                VectorXd vec_in_list;
                get_vector_from_bed_matrix(MACOJO::final_commonSNP_index[SNP_list[j]], vec_in_list);
                r_temp_vec(j) = (X_temp_vec.transpose() * vec_in_list).value() / (indi_num-1);
            }       
        }
    }
}


void Cohort::calc_conditional_effects() 
{   
    MatrixXd temp1, temp2;
    ArrayXd temp3;

    temp1 = sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5));
    temp2.noalias() = R_inv_pre * temp1;
    temp3 = r * temp2;
    conditional_beta = sumstat_screened.col(0) - temp3 / sqrt(sumstat_screened.col(5));
}


bool Cohort::calc_joint_effects(const ArrayXXd &s, bool flag, double iter_colinear_threshold) 
{      
    if (flag && ((abs(R_inv_post.minCoeff()) > iter_colinear_threshold) || (abs(R_inv_post.maxCoeff()) > iter_colinear_threshold)))
        return true;

    VectorXd temp1 = sqrt(s.col(5)) * s.col(0);
    ArrayXd temp2 = R_inv_post * temp1;
    beta = temp2 / sqrt(s.col(5));

    double sigma_J_squared = Vp - (s.col(0) * s.col(5) * beta).sum();
                
    beta_var = sigma_J_squared * R_inv_post.diagonal().array() / s.col(6) ;

    if (flag && beta_var.minCoeff() <= 0)
        return true;

    double Neff = median(s.col(4));
    int M = s.rows();
    R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    return false;
}


void Cohort::calc_R_inv(bool if_fast_inv) {
    if (if_fast_inv) {
        VectorXd temp_vector = R_inv_pre * r_temp_vec;
        double temp_element = 1 / (1 - r_temp_vec.transpose() * temp_vector);
    
        int dim = R_inv_pre.rows();
        R_inv_post.setZero(dim+1, dim+1);
    
        R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
        R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
        R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
        R_inv_post(dim, dim) = temp_element;
    } else {
        int dim = R_pre.rows();
        R_post.setZero(dim+1, dim+1);
    
        R_post.block(0, 0, dim, dim) = R_pre;
        R_post.block(0, dim, dim, 1) = r_temp_vec;
        R_post.block(dim, 0, 1, dim) = r_temp_vec.transpose();
        R_post(dim, dim) = 1;

        R_inv_post.noalias() = R_post.ldlt().solve(MatrixXd::Identity(dim+1, dim+1));
    }
}


void Cohort::calc_R_inv_from_SNP_list(const vector<int> &SNP_list, int window_size, bool if_LD_mode) 
{   
    int total_num = SNP_list.size();
    R_post = MatrixXd::Identity(total_num, total_num);

    if (if_LD_mode) {
        for (int i = 0; i < total_num; i++) {
            for (int j = 0; j < i; j++) {
                if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[i]]] - 
                        MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]]) <= window_size) {
                    R_post(i, j) = LD_matrix(MACOJO::final_commonSNP_index[SNP_list[i]], MACOJO::final_commonSNP_index[SNP_list[j]]);
                    R_post(j, i) = R_post(i, j);
                }
            }
        }
    } else {
        vector<VectorXd> X_new_model;
        
        for (int index : SNP_list) {
            get_vector_from_bed_matrix(MACOJO::final_commonSNP_index[index], X_temp_vec);
            X_new_model.push_back(X_temp_vec);
        }

        for (int i = 0; i < total_num; i++) {
            for (int j = 0; j < i; j++) {
                if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[i]]] - 
                        MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]]) <= window_size) {
                    R_post(i, j) = (X_new_model[i].transpose() * X_new_model[j]).value() / (indi_num-1);
                    R_post(j, i) = R_post(i, j);
                }
            }
        }
    } 

    // R_inv_post.noalias() = R_post.inverse();
    R_inv_post.noalias() = R_post.ldlt().solve(MatrixXd::Identity(total_num, total_num));
}
    

void MACOJO::accept_SNP_as_candidate(int candidate_index) 
{   
    candidate_SNP.push_back(candidate_index);

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];
        c.calc_inner_product_with_SNP_list(screened_SNP, candidate_index, window_size, if_LD_mode);
        append_column(c.r, c.r_temp_vec);
    }
}


void MACOJO::remove_new_colinear_SNP(int candidate_index) 
{   
    for (int i = screened_SNP.size()-1, last_col = candidate_SNP.size()-1; i >= 0; i--) {
        
        bool erase_flag = false;
        for (int n : current_calculation_list) {
            if (abs(cohorts[n].r(i, last_col)) >= colinear_threshold_sqrt) {
                erase_flag = true; 
                break;
            }
        }

        if (erase_flag) {
            for (int n : current_calculation_list) {
                remove_row(cohorts[n].sumstat_screened, i);
                remove_row(cohorts[n].r, i);
            }
            
            excluded_SNP.push_back(screened_SNP[i]);
            screened_SNP.erase(screened_SNP.begin()+i);
        }
    }

    auto iter = find(excluded_SNP.begin(), excluded_SNP.end(), candidate_index);
    if (iter != excluded_SNP.end())
        excluded_SNP.erase(iter);
}   


void MACOJO::initialize_MDISA_from_MACOJO(Cohort &c) 
{   
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {
        c.calc_inner_product_with_SNP_list(candidate_SNP, excluded_SNP[i], window_size, if_LD_mode);
        
        if (c.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) {
            append_row(c.sumstat_screened, c.sumstat.row(excluded_SNP[i]));
            append_row(c.r, c.r_temp_vec.transpose());
            screened_SNP.push_back(excluded_SNP[i]); 
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }

    for (int i = backward_removed_SNP.size()-1; i >= 0; i--) {
        c.calc_inner_product_with_SNP_list(candidate_SNP, backward_removed_SNP[i], window_size, if_LD_mode);
        
        if (c.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) {
            append_row(c.sumstat_screened, c.sumstat.row(backward_removed_SNP[i])); 
            append_row(c.r, c.r_temp_vec.transpose());
            screened_SNP.push_back(backward_removed_SNP[i]);
        } else {
            excluded_SNP.push_back(backward_removed_SNP[i]); 
        }
    }

    vector<int>().swap(backward_removed_SNP);
}


void MACOJO::adjust_SNP_according_to_backward_selection() 
{   
    for (int i = candidate_SNP.size()-1; i >= fixed_candidate_SNP_num; i--) {

        if (sumstat_new_model_joint(i, 3) > threshold) {
            for (int n : current_calculation_list) {
                auto &c = cohorts[n];
                remove_row(c.sumstat_candidate, i);
                remove_column(c.r, i);
            }
                
            LOGGER.w(1, "Previous candidate SNP removed", final_commonSNP[candidate_SNP[i]]);
            backward_removed_SNP.push_back(candidate_SNP[i]);
            candidate_SNP.erase(candidate_SNP.begin()+i);
        }
    }
    
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {

        bool append_flag = true;
        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            // check if the excluded SNP is colinear with any remaining candidate SNP
            c.calc_inner_product_with_SNP_list(candidate_SNP, excluded_SNP[i], window_size, if_LD_mode);
            if (c.r_temp_vec.cwiseAbs().maxCoeff() >= colinear_threshold_sqrt) {
                append_flag = false; 
                break;
            }
        }
    
        if (append_flag) {
            for (int n : current_calculation_list) {
                auto &c = cohorts[n];
                append_row(c.sumstat_screened, c.sumstat.row(excluded_SNP[i]));
                append_row(c.r, c.r_temp_vec.transpose());
            }
            
            screened_SNP.push_back(excluded_SNP[i]); 
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }
}


void MACOJO::save_temp_model() 
{   
    for (int n : current_calculation_list) {
        auto &c = cohorts[n];

        c.R_inv_pre = c.R_inv_post;

        // only used for if_fast_inv == false, saved anyway to avoid complexity
        c.R_pre = c.R_post;

        c.previous_R2 = c.R2;
        c.output_b = c.beta;
        c.output_se2 = c.beta_var;
    }
}


void MACOJO::inverse_var_meta(bool if_conditional) 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    sumstat_merge.setZero(screened_SNP.size(), 4);

    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        sumstat_merge.col(0) = if_conditional ? c.conditional_beta : c.sumstat_screened.col(0);
        sumstat_merge.col(1) = c.sumstat_screened.col(1);
    } else {
        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]];
            sumstat_merge.col(0) += (if_conditional ? c.conditional_beta : c.sumstat_screened.col(0)) / c.sumstat_screened.col(1);
            sumstat_merge.col(1) += 1 / c.sumstat_screened.col(1);
        }

        sumstat_merge.col(1) = 1 / sumstat_merge.col(1);
        sumstat_merge.col(0) = sumstat_merge.col(0) * sumstat_merge.col(1);
    }

    sumstat_merge.col(2) = abs(sumstat_merge.col(0) / sqrt(sumstat_merge.col(1)));
    sumstat_merge.col(3) = erfc(sumstat_merge.col(2) / sqrt(2));
}


void MACOJO::inverse_var_meta_joint() 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    sumstat_new_model_joint.setZero(candidate_SNP.size(), 4);

    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        sumstat_new_model_joint.col(0) = c.beta;
        sumstat_new_model_joint.col(1) = c.beta_var;
    } else {
        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]];
            sumstat_new_model_joint.col(0) += c.beta / c.beta_var;
            sumstat_new_model_joint.col(1) += 1 / c.beta_var;
        }

        sumstat_new_model_joint.col(1) = 1 / sumstat_new_model_joint.col(1);
        sumstat_new_model_joint.col(0) = sumstat_new_model_joint.col(0) * sumstat_new_model_joint.col(1);
    }

    sumstat_new_model_joint.col(2) = abs(sumstat_new_model_joint.col(0) / sqrt(sumstat_new_model_joint.col(1)));
    sumstat_new_model_joint.col(3) = erfc(sumstat_new_model_joint.col(2) / sqrt(2));
}
