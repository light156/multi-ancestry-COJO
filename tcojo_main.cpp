#include "tcojo.h"


void TCOJO::main_loop(string savename) 
{   
    inverse_var_meta(c1.sumstat_screened.col(0), c2.sumstat_screened.col(0), \
        c1.sumstat_screened.col(1), c2.sumstat_screened.col(1), sumstat_merge);
    
    sumstat_merge.col(2).maxCoeff(&max_SNP_index);
    if (sumstat_merge(max_SNP_index, 3) > threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    candidate_SNP.push_back(max_SNP_index);
    screened_SNP.resize(final_commonSNP.size());
    iota(screened_SNP.begin(), screened_SNP.end(), 0);

    LOGGER.i(0, "First SNP", final_commonSNP[max_SNP_index]);
    LOGGER << "--------------------------------" << endl;

    initialize_matrices(c1);
    initialize_matrices(c2);
    remove_new_colinear_SNP();

    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " "
            << excluded_SNP.size() << " " << backward_removed_SNP.size() << endl;

        // calculate conditional effects
        c1.calc_conditional_effects();
        c2.calc_conditional_effects();
        inverse_var_meta(c1.conditional_beta, c2.conditional_beta, 
            c1.sumstat_screened.col(1), c2.sumstat_screened.col(1), sumstat_merge);
        
        while (true) {
            // select maximal SNP
            double max_Zabs = sumstat_merge.col(2).maxCoeff(&max_SNP_index);
            
            if (max_Zabs < 0 || sumstat_merge(max_SNP_index, 3) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }
            
            string max_SNP_name = final_commonSNP[screened_SNP[max_SNP_index]];

            // calculate joint effects Cohort 1
            append_row(c1.sumstat_candidate, c1.sumstat_screened.row(max_SNP_index));

            calc_inner_product_with_candidate(c1, screened_SNP[max_SNP_index]);
            c1.calc_R_inv(if_fast_inv);

            NA_flag = c1.calc_joint_effects(c1.sumstat_candidate, true, iter_colinear_threshold);
            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            // calculate joint effects Cohort 2
            append_row(c2.sumstat_candidate, c2.sumstat_screened.row(max_SNP_index));

            calc_inner_product_with_candidate(c2, screened_SNP[max_SNP_index]);
            c2.calc_R_inv(if_fast_inv);

            NA_flag = c2.calc_joint_effects(c2.sumstat_candidate, true, iter_colinear_threshold);
            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                remove_row(c2.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            if ((c1.R2 < (1+R2_incremental_threshold) * c1.previous_R2) || 
                (c2.R2 < (1+R2_incremental_threshold) * c2.previous_R2)) {
                // LOGGER.w(1, "R2 increment unsatisfactory", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                remove_row(c2.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            // temporarily include new candidate SNP
            candidate_SNP.push_back(screened_SNP[max_SNP_index]);

            c1.X_candidate.push_back(c1.X_temp_vec);
            calc_inner_product_with_screened(c1, screened_SNP[max_SNP_index]);
            append_column(c1.r, c1.r_temp_vec);
            
            c2.X_candidate.push_back(c2.X_temp_vec);
            calc_inner_product_with_screened(c2, screened_SNP[max_SNP_index]);
            append_column(c2.r, c2.r_temp_vec);
            
            // check joint model threshold
            inverse_var_meta(c1.beta, c2.beta, c1.beta_var, c2.beta_var, sumstat_new_model_joint);

            if (sumstat_new_model_joint.col(3).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed", max_SNP_name);

                int M = candidate_SNP.size();
                // LOGGER << "Added R vector cohort 1: " << c1.r.row(max_SNP_index) << endl;
                // LOGGER << "Added R vector cohort 2: " << c2.r.row(max_SNP_index) << endl;
                remove_new_colinear_SNP();
                LOGGER << "Added diagonal value cohort 1: " << c1.R_inv_post(M-1, M-1) << endl;
                LOGGER << "Added diagonal value cohort 2: " << c2.R_inv_post(M-1, M-1) << endl;
                LOGGER << "Joint b: " << sumstat_new_model_joint(M-1, 0) << endl;
                LOGGER << "Joint se: " << sqrt(sumstat_new_model_joint(M-1, 1)) << endl;
                LOGGER << "Joint p-value: " << scientific << sumstat_new_model_joint(M-1, 3) << endl;
                LOGGER << "Adjusted R2 for cohort 1: " << fixed << c1.R2 << endl;
                LOGGER << "Adjusted R2 for cohort 2: " << c2.R2 << endl;
                break; 
            }

            // backward selection
            LOGGER << "Backward selection" << endl;

            initialize_backward_selection(c1, sumstat_new_model_joint.col(3));
            initialize_backward_selection(c2, sumstat_new_model_joint.col(3));

            c1.calc_joint_effects(c1.sumstat_backward_new_model, false);
            c2.calc_joint_effects(c2.sumstat_backward_new_model, false);

            LOGGER << "Adjusted R2 after SNP elimination for cohort 1: " << fixed << c1.R2 << endl;
            LOGGER << "Adjusted R2 after SNP elimination for cohort 2: " << c2.R2 << endl;

            if ((c1.R2 < (1+R2_incremental_threshold_backwards) * c1.previous_R2) || 
                (c2.R2 < (1+R2_incremental_threshold_backwards) * c2.previous_R2)) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", max_SNP_name);
                candidate_SNP.pop_back();

                remove_row(c1.sumstat_candidate);
                c1.X_candidate.pop_back();
                remove_column(c1.r);

                remove_row(c2.sumstat_candidate);
                c2.X_candidate.pop_back();
                remove_column(c2.r);
                
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", max_SNP_name);
            adjust_SNP_according_to_backward_selection(sumstat_new_model_joint.col(3));
            remove_new_colinear_SNP();
            break;
        }

        // save template model for output
        if (!loop_break_indicator) {
            c1.save_temp_model();
            c2.save_temp_model();
        }

        LOGGER << "iter " << ++iter_num << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }

    inverse_var_meta(c1.output_b, c2.output_b, c1.output_se2, c2.output_se2, sumstat_merge);
    save_results_main_loop(savename + ".jma.cojo");

    // backup SNPs for Cohort 1
    MDISA_fixed_candidate_SNP_num = candidate_SNP.size();
    
    vector<int> candidate_SNP_backup = candidate_SNP;
    vector<int> screened_SNP_backup = screened_SNP;
    vector<int> excluded_SNP_backup = excluded_SNP;
    vector<int> backward_removed_SNP_backup = backward_removed_SNP;

    // MDISA Cohort 1
    LOGGER.i(0, "Cohort 1", "MDISA");
    initialize_MDISA(c1);
    MDISA(c1);
    save_results_DISA(c1, savename + ".MDISA.cohort1.jma.cojo");

    // backup SNPs for Cohort 2
    candidate_SNP = candidate_SNP_backup;
    screened_SNP = screened_SNP_backup;
    excluded_SNP = excluded_SNP_backup;
    backward_removed_SNP = backward_removed_SNP_backup;

    // MDISA Cohort 2
    LOGGER.i(0, "Cohort 2", "MDISA");
    initialize_MDISA(c2);
    MDISA(c2);
    save_results_DISA(c2, savename + ".MDISA.cohort2.jma.cojo");

    vector<int>().swap(candidate_SNP_backup);
    vector<int>().swap(screened_SNP_backup);
    vector<int>().swap(excluded_SNP_backup);
    vector<int>().swap(backward_removed_SNP_backup);
}


void TCOJO::MDISA(Cohort &c) 
{   
    bool cohort1_only = (addressof(c) == addressof(c1));
    ArrayXd Zabs_temp, p_temp, ZabsJ, pJ;

    if (candidate_SNP.size() == 0) {

        Zabs_temp = abs(c.sumstat_screened.col(0) / c.sumstat_screened.col(1));
        p_temp = erfc(Zabs_temp / sqrt(2));

        Zabs_temp.maxCoeff(&max_SNP_index);
        if (p_temp(max_SNP_index) > threshold)
            LOGGER.e(0, "Input data has no significant SNPs.");
        
        candidate_SNP.push_back(max_SNP_index);
        screened_SNP.resize(final_commonSNP.size());
        iota(screened_SNP.begin(), screened_SNP.end(), 0);

        LOGGER.i(0, "First SNP", final_commonSNP[max_SNP_index]);
        LOGGER << "--------------------------------" << endl;

        // impossible to happen, just in case
        vector<int>().swap(excluded_SNP);
        vector<int>().swap(backward_removed_SNP);
        
        initialize_matrices(c);
        remove_new_colinear_SNP(cohort1_only, !cohort1_only);
    } 

    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " "
            << excluded_SNP.size() << " " << backward_removed_SNP.size() << endl;

        // calculate conditional effects
        c.calc_conditional_effects();
        Zabs_temp = abs(c.conditional_beta / sqrt(c.sumstat_screened.col(1)));
        p_temp = erfc(Zabs_temp / sqrt(2));
    
        while (true) {
            // select maximal SNP
            double max_Zabs = Zabs_temp.maxCoeff(&max_SNP_index);

            if (max_Zabs < 0 || p_temp(max_SNP_index) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }
            
            string max_SNP_name = final_commonSNP[screened_SNP[max_SNP_index]];

            // calculate joint effects
            append_row(c.sumstat_candidate, c.sumstat_screened.row(max_SNP_index));

            calc_inner_product_with_candidate(c, screened_SNP[max_SNP_index]);
            c.calc_R_inv(if_fast_inv);

            NA_flag = c.calc_joint_effects(c.sumstat_candidate, true, iter_colinear_threshold);
            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c.sumstat_candidate);
                Zabs_temp(max_SNP_index) = -1;
                continue;
            }

            if (c.R2 < (1+R2_incremental_threshold) * c.previous_R2) {
                // LOGGER.w(1, "R2 increment unsatisfactory", max_SNP_name);
                remove_row(c.sumstat_candidate);
                Zabs_temp(max_SNP_index) = -1;
                continue;
            }

            // temporarily include new candidate SNP
            candidate_SNP.push_back(screened_SNP[max_SNP_index]);

            c.X_candidate.push_back(c.X_temp_vec);
            calc_inner_product_with_screened(c, screened_SNP[max_SNP_index]);
            append_column(c.r, c.r_temp_vec);
            
            // check joint model threshold
            ZabsJ = abs(c.beta / sqrt(c.beta_var));
            pJ = erfc(ZabsJ / sqrt(2));

            if (pJ.bottomRows(candidate_SNP.size()-MDISA_fixed_candidate_SNP_num).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed", max_SNP_name);

                int M = candidate_SNP.size();
                // LOGGER << "Added R vector: " << c.r.row(max_SNP_index) << endl;
                remove_new_colinear_SNP(cohort1_only, !cohort1_only);
                LOGGER << "Added diagonal value: " << c.R_inv_post(M-1, M-1) << endl;
                LOGGER << "Joint b: " << c.beta(M-1) << endl;
                LOGGER << "Joint se: " << sqrt(c.beta_var(M-1)) << endl;
                LOGGER << "Joint p-value: " << scientific << pJ(M-1) << endl;
                LOGGER << "Adjusted R2: " << fixed << c.R2 << endl;
                break; 
            }

            // backward selection
            LOGGER << "Backward selection" << endl;

            initialize_backward_selection(c, pJ);
            c.calc_joint_effects(c.sumstat_backward_new_model, false);

            LOGGER << "Adjusted R2 after SNP elimination: " << fixed << c.R2 << endl;

            if (c.R2 < (1+R2_incremental_threshold_backwards) * c.previous_R2) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", max_SNP_name);
                candidate_SNP.pop_back();

                remove_row(c.sumstat_candidate);
                c.X_candidate.pop_back();
                remove_column(c.r);

                Zabs_temp(max_SNP_index) = -1;
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", max_SNP_name);
            adjust_SNP_according_to_backward_selection(pJ, cohort1_only, !cohort1_only);
            remove_new_colinear_SNP(cohort1_only, !cohort1_only);
            break;
        }

        // save template model for output
        if (!loop_break_indicator)
            c.save_temp_model();

        LOGGER << "iter " << ++iter_num << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }
}