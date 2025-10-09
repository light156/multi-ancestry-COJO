#include "macojo.h"


void MACOJO::accept_SNP_as_candidate(int screened_index) 
{   
    int candidate_index = screened_SNP[screened_index];
    candidate_SNP.push_back(candidate_index);
    screened_SNP.erase(screened_SNP.begin() + screened_index);

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];
        
        // save current model
        if (candidate_SNP.size() == 1) {
            c.output_b = c.sumstat_candidate.col(0);
            c.output_se2 = c.sumstat_candidate.col(1);
            if (if_gcta_COJO) {
                c.R_pre(0,0) = c.sumstat_candidate(0,6);
                c.R_inv_pre(0,0) = 1.0 / c.sumstat_candidate(0,6);
            }
        } else {
            c.output_b = c.beta;
            c.output_se2 = c.beta_var;
            c.R_inv_pre = c.R_inv_post;
            c.R_pre = c.R_post; // only used for if_fast_inv == false, saved anyway to avoid complexity
            c.previous_R2 = c.R2; // only used for if_gcta_COJO == false, saved anyway to avoid complexity
        }

        remove_row(c.sumstat_screened, screened_index);
        remove_row(c.r, screened_index);
        c.calc_inner_product_with_SNP_list(screened_SNP, c.sumstat_screened, candidate_index);
        append_column(c.r, c.r_temp_vec);
    }
}


void MACOJO::adjust_SNP_according_to_backward_selection() 
{   
    for (int i = sumstat_new_model_joint.rows()-1; i >= fixed_candidate_SNP_num; i--) {

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
}


void MACOJO::initialize_MDISA_from_MACOJO(Cohort &c) 
{   
    for (int i = backward_removed_SNP.size()-1; i >= 0; i--) {
        c.calc_inner_product_with_SNP_list(candidate_SNP, c.sumstat_candidate, backward_removed_SNP[i]);
        
        if (c.r_temp_vec.cwiseAbs().maxCoeff() < sqrt(colinear_threshold)) {
            append_row(c.sumstat_screened, c.sumstat.row(backward_removed_SNP[i]));
            append_row(c.r, c.r_temp_vec.transpose());
            screened_SNP.push_back(backward_removed_SNP[i]); 
            backward_removed_SNP.erase(backward_removed_SNP.begin()+i);
        }
    }
}


void MACOJO::inverse_var_meta_init() 
{
    // merge: 2:Zabs, 3:p
    sumstat_merge.setZero(screened_SNP.size(), 4);

    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        sumstat_merge.col(2) = abs(c.sumstat_screened.col(0) / sqrt(c.sumstat_screened.col(1)));
        sumstat_merge.col(3) = erfc(sumstat_merge.col(2) / sqrt(2));
    } else {
        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]];
            sumstat_merge.col(0) += c.sumstat_screened.col(0) / c.sumstat_screened.col(1);
            sumstat_merge.col(1) += 1 / c.sumstat_screened.col(1);
        }
        sumstat_merge.col(2) = abs(sumstat_merge.col(0)) / sqrt(sumstat_merge.col(1));
        sumstat_merge.col(3) = erfc(sumstat_merge.col(2) / sqrt(2));
    }
}


void MACOJO::inverse_var_meta_conditional() 
{
    // merge: 2:Zabs, 3:p
    sumstat_merge.setZero(screened_SNP.size(), 4);

    if (current_calculation_list.size() == 1) {
        auto &c = cohorts[current_calculation_list[0]];
        if (if_gcta_COJO)
            sumstat_merge.col(2) = abs(c.conditional_beta) * sqrt(c.sumstat_screened.col(6)) / sqrt(c.Vp);
        else
            sumstat_merge.col(2) = abs(c.conditional_beta) / sqrt(c.sumstat_screened.col(1));

        sumstat_merge.col(3) = erfc(sumstat_merge.col(2) / sqrt(2));
    } else {
        for (int n : current_calculation_list) {
            auto &c = cohorts[current_calculation_list[n]];
            if (if_gcta_COJO) {
                sumstat_merge.col(0) += c.conditional_beta * c.sumstat_screened.col(6) / c.Vp;
                sumstat_merge.col(1) += c.sumstat_screened.col(6) / c.Vp;
            }
            else {
                sumstat_merge.col(0) += c.conditional_beta / c.sumstat_screened.col(1);
                sumstat_merge.col(1) += 1 / c.sumstat_screened.col(1);
            }
        }   
        sumstat_merge.col(2) = abs(sumstat_merge.col(0)) / sqrt(sumstat_merge.col(1));
        sumstat_merge.col(3) = erfc(sumstat_merge.col(2) / sqrt(2));
    }
}


void MACOJO::inverse_var_meta_joint() 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    sumstat_new_model_joint.setZero(candidate_SNP.size()+1, 4);

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

    sumstat_new_model_joint.col(2) = abs(sumstat_new_model_joint.col(0)) / sqrt(sumstat_new_model_joint.col(1));
    sumstat_new_model_joint.col(3) = erfc(sumstat_new_model_joint.col(2) / sqrt(2));
}
