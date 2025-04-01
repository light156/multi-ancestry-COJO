#include "tcojo.h"


void TCOJO::calc_inner_product_with_candidate(Cohort &c, int single_index)
{   
    c.get_vector_from_bed_matrix(final_commonSNP_index[single_index], c.X_temp_vec);
    int X_temp_vec_pos_ref = SNP_pos_ref[final_commonSNP_index[single_index]];

    c.r_temp_vec.setZero(candidate_SNP.size());
    
    #pragma omp parallel for
    for (int j = 0; j < candidate_SNP.size(); j++) {
        if (abs(SNP_pos_ref[final_commonSNP_index[candidate_SNP[j]]] - X_temp_vec_pos_ref) <= window_size) {
            c.r_temp_vec(j) = (c.X_temp_vec.transpose() * c.X_candidate[j]).value() / (c.indi_num-1);
        }       
    }
}


// must be called after get_vector_from_bed_matrix()
void TCOJO::calc_inner_product_with_screened(Cohort &c, int single_index)
{   
    int X_temp_vec_pos_ref = SNP_pos_ref[final_commonSNP_index[single_index]];
    c.r_temp_vec.setZero(screened_SNP.size());

    VectorXd vec_in_list;

    #pragma omp parallel for 
    for (int j = 0; j < screened_SNP.size(); j++) {
        if (abs(SNP_pos_ref[final_commonSNP_index[screened_SNP[j]]] - X_temp_vec_pos_ref) <= window_size) {
            c.get_vector_from_bed_matrix(final_commonSNP_index[screened_SNP[j]], vec_in_list);
            c.r_temp_vec(j) = (c.X_temp_vec.transpose() * vec_in_list).value() / (c.indi_num-1);
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
        R_inv_post.resize(dim+1, dim+1);
    
        R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
        R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
        R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
        R_inv_post(dim, dim) = temp_element;
    } else {
        int dim = R_pre.rows();
        R_post.resize(dim+1, dim+1);
    
        R_post.block(0, 0, dim, dim) = R_pre;
        R_post.block(0, dim, dim, 1) = r_temp_vec;
        R_post.block(dim, 0, 1, dim) = r_temp_vec.transpose();
        R_post(dim, dim) = 1;

        R_inv_post.noalias() = R_post.ldlt().solve(MatrixXd::Identity(dim+1, dim+1));
    }
}


void Cohort::save_temp_model() 
{   
    R_inv_pre = R_inv_post;

    // only used for if_fast_inv == false, saved anyway to avoid complexity
    R_pre = R_post;

    previous_R2 = R2;
    output_b = beta;
    output_se2 = beta_var;
}


void TCOJO::inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
    const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge) 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    merge.resize(b_cohort1.size(), 4);

    merge.col(0) = (b_cohort1 * se2_cohort2 + b_cohort2 * se2_cohort1) / (se2_cohort1 + se2_cohort2);
    merge.col(1) = se2_cohort1 * se2_cohort2 / (se2_cohort1 + se2_cohort2);
    merge.col(2) = abs(merge.col(0) / sqrt(merge.col(1)));
    merge.col(3) = erfc(merge.col(2) / sqrt(2));
}


void TCOJO::initialize_matrices(Cohort &c) 
{   
    c.sumstat_candidate = c.sumstat_screened.row(max_SNP_index);
    
    c.get_vector_from_bed_matrix(final_commonSNP_index[max_SNP_index], c.X_temp_vec);
    c.X_candidate.push_back(c.X_temp_vec);

    calc_inner_product_with_screened(c, max_SNP_index);
    c.r = c.r_temp_vec;
    
    c.R_inv_pre = MatrixXd::Identity(1,1);
    c.R_pre = MatrixXd::Identity(1,1);
    
    c.previous_R2 = 0.0;
    c.output_b = c.sumstat_screened.col(0);
    c.output_se2 = c.sumstat_screened.col(1);
}


void TCOJO::initialize_MDISA(Cohort &c) 
{   
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {
        calc_inner_product_with_candidate(c, excluded_SNP[i]); 
        
        if (c.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) {
            append_row(c.sumstat_screened, c.sumstat.row(excluded_SNP[i]));
            append_row(c.r, c.r_temp_vec.transpose());
            screened_SNP.push_back(excluded_SNP[i]); 
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }

    for (int i = backward_removed_SNP.size()-1; i >= 0; i--) {
        calc_inner_product_with_candidate(c, backward_removed_SNP[i]);
        
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


void TCOJO::initialize_backward_selection(Cohort &c, const ArrayXd &pJ)
{   
    c.sumstat_backward_new_model = c.sumstat_candidate;
    c.X_backward_new_model = c.X_candidate;

    for (int i = candidate_SNP.size()-1; i >= MDISA_fixed_candidate_SNP_num; i--) {
        if (pJ(i) > threshold) {
            remove_row(c.sumstat_backward_new_model, i);
            c.X_backward_new_model.erase(c.X_backward_new_model.begin()+i);
        }
    }

    int total_num = c.X_backward_new_model.size();
    c.R_post = MatrixXd::Identity(total_num, total_num);

    vector<int> SNP_pos_backward;
    for (int i = 0; i < candidate_SNP.size(); i++) {
        if (i < MDISA_fixed_candidate_SNP_num || pJ(i) <= threshold)
            SNP_pos_backward.push_back(SNP_pos_ref[final_commonSNP_index[candidate_SNP[i]]]);
    }

    #pragma omp parallel for 
    for (int i = 1; i < total_num; i++)
        for (int j = 0; j < i; j++)
            if (abs(SNP_pos_backward[i] - SNP_pos_backward[j]) <= window_size) {
                c.R_post(i, j) = (c.X_backward_new_model[i].transpose() * c.X_backward_new_model[j]).value() / (c.indi_num-1);
                c.R_post(j, i) = c.R_post(i, j);
            } 
    
    vector<int>().swap(SNP_pos_backward);
    
    if (if_fast_inv)
        c.R_inv_post.noalias() = c.R_post.inverse();
    else
        c.R_inv_post.noalias() = c.R_post.ldlt().solve(MatrixXd::Identity(total_num, total_num));
}


void TCOJO::remove_new_colinear_SNP(bool cohort1_only, bool cohort2_only) 
{   
    for (int i = screened_SNP.size()-1, last_col = candidate_SNP.size()-1; i >= 0; i--) {
        if ((!cohort2_only && abs(c1.r(i, last_col)) >= colinear_threshold_sqrt) || 
            (!cohort1_only && abs(c2.r(i, last_col)) >= colinear_threshold_sqrt)) {
            if (!cohort2_only) {
                remove_row(c1.sumstat_screened, i);
                remove_row(c1.r, i);
            }

            if (!cohort1_only) {
                remove_row(c2.sumstat_screened, i);
                remove_row(c2.r, i);
            }
            
            if (i != max_SNP_index)
                excluded_SNP.push_back(screened_SNP[i]);
            screened_SNP.erase(screened_SNP.begin()+i);
        }
    }
}   


void TCOJO::adjust_SNP_according_to_backward_selection(const ArrayXd &pJ, bool cohort1_only, bool cohort2_only) 
{   
    if (!cohort2_only) {
        remove_row(c1.sumstat_screened, max_SNP_index);
        remove_row(c1.r, max_SNP_index);
        c1.sumstat_backward_new_model.swap(c1.sumstat_candidate);
        c1.X_backward_new_model.swap(c1.X_candidate);
    }

    if (!cohort1_only) {
        remove_row(c2.sumstat_screened, max_SNP_index);
        remove_row(c2.r, max_SNP_index);
        c2.sumstat_backward_new_model.swap(c2.sumstat_candidate);
        c2.X_backward_new_model.swap(c2.X_candidate);
    }

    screened_SNP.erase(screened_SNP.begin()+max_SNP_index);

    for (int i = candidate_SNP.size()-1; i >= MDISA_fixed_candidate_SNP_num; i--) {
        if (pJ(i) > threshold) {
            if (!cohort2_only) 
                remove_column(c1.r, i);
            
            if (!cohort1_only)
                remove_column(c2.r, i);

            LOGGER.w(1, "Previous candidate SNP removed", final_commonSNP[candidate_SNP[i]]);
            backward_removed_SNP.push_back(candidate_SNP[i]);
            candidate_SNP.erase(candidate_SNP.begin()+i);
        }
    }
    
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {
        if (!cohort2_only)
            calc_inner_product_with_candidate(c1, excluded_SNP[i]);
    
        if (!cohort1_only)
            calc_inner_product_with_candidate(c2, excluded_SNP[i]);
        
        if ((cohort2_only || c1.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) && 
            (cohort1_only || c2.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt)) {
            if (!cohort2_only) {
                append_row(c1.sumstat_screened, c1.sumstat.row(excluded_SNP[i]));
                append_row(c1.r, c1.r_temp_vec.transpose());
            }

            if (!cohort1_only) {
                append_row(c2.sumstat_screened, c2.sumstat.row(excluded_SNP[i])); 
                append_row(c2.r, c2.r_temp_vec.transpose());
            }
            
            screened_SNP.push_back(excluded_SNP[i]); 
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }
}
