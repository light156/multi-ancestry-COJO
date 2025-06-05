#include "macojo.h"


void MACOJO::entry_function(string savename)
{   
    LOGGER.i(0, "Calculation started!");

    for (int n = 0; n < cohorts.size(); n++)
        current_calculation_list.push_back(n);

    initialize_main_loop();

    if (if_cojo_joint && fixed_candidate_SNP_num > 0) {
        // output results with fixed candidate SNPs and exit
        output_results_to_file(savename+".jma.cojo.fixedSNP");
        LOGGER.i(0, "Calculation finished with fixed candidate SNPs, program exit!");
        return;
    }

    main_loop();

    // handle case for single cohort
    if (cohorts.size() == 1) {
        output_results_to_file(savename+".MDISA.jma.cojo");
        return;
    }

    // save results for MACOJO
    output_results_to_file(savename+".jma.cojo");

    if (if_MDISA) {
        // Start MDISA analysis, backup SNPs for MDISA
        fixed_candidate_SNP_num = candidate_SNP.size();
        vector<int> candidate_SNP_backup(candidate_SNP);
        vector<int> screened_SNP_backup(screened_SNP);
        vector<int> excluded_SNP_backup(excluded_SNP);
        vector<int> backward_removed_SNP_backup(backward_removed_SNP);

        current_calculation_list.resize(1);

        for (int temp_index = 0; temp_index < cohorts.size(); temp_index++) {
            LOGGER << "MDISA for Cohort " << temp_index+1 << endl;
            current_calculation_list[0] = temp_index;

            if (temp_index>0) {                 
                candidate_SNP = candidate_SNP_backup;
                screened_SNP = screened_SNP_backup;
                excluded_SNP = excluded_SNP_backup;
                backward_removed_SNP = backward_removed_SNP_backup;
            }

            initialize_MDISA_from_MACOJO(cohorts[temp_index]);
            main_loop();
            output_results_to_file(savename+".MDISA.cohort" + to_string(temp_index+1) + ".jma.cojo");
        }
    }

    LOGGER.i(0, "All calculations finished!");
}


void MACOJO::initialize_main_loop() 
{   
    screened_SNP.resize(final_commonSNP.size());
    for (int i = 0; i < final_commonSNP.size(); i++)
        screened_SNP[i] = i;

    // deal with cases where user provides fixed candidate SNPs
    for (const auto &snp_name : fixed_candidate_SNP) {
        int single_index = distance(final_commonSNP.begin(), find(final_commonSNP.begin(), final_commonSNP.end(), snp_name));

        if (single_index < final_commonSNP.size()) {
            for (int n : current_calculation_list) 
                append_row(cohorts[n].sumstat_candidate, cohorts[n].sumstat.row(single_index));
            
            if (!if_cojo_joint) {
                accept_SNP_as_candidate(single_index);
                remove_new_colinear_SNP(single_index);
            } else 
                candidate_SNP.push_back(single_index);
        }
    }

    fixed_candidate_SNP_num = candidate_SNP.size();

    if (fixed_candidate_SNP_num > 0) {
        LOGGER << fixed_candidate_SNP_num << " fixed candidate SNPs provided by the user" << endl;

        // adjust for case n==1
        if (fixed_candidate_SNP_num == 1) {
            for (int n : current_calculation_list) {
                auto &c = cohorts[n];
                c.output_b = c.sumstat_candidate.col(0);
                c.output_se2 = c.sumstat_candidate.col(1);
            }
            return;
        }

        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            c.calc_R_inv_from_SNP_list(candidate_SNP, window_size);
            c.calc_joint_effects(c.sumstat_candidate, false);
        }

        save_temp_model();
    }
}


void MACOJO::main_loop() 
{   
    int max_screened_index, true_index;

    if (candidate_SNP.size() == 0) {
        LOGGER.i(0, "No candidate SNPs, using the most significant SNP as the first candidate");

        inverse_var_meta(false);

        sumstat_merge.col(2).maxCoeff(&max_screened_index);
        if (sumstat_merge(max_screened_index, 3) > threshold)
            LOGGER.e(0, "Input data has no significant SNPs");

        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            c.sumstat_candidate = c.sumstat_screened.row(max_screened_index);
            c.output_b = c.sumstat_candidate.col(0);
            c.output_se2 = c.sumstat_candidate.col(1);
        }
        
        true_index = screened_SNP[max_screened_index];
        accept_SNP_as_candidate(true_index);
        remove_new_colinear_SNP(true_index);

        LOGGER.i(0, "First SNP", final_commonSNP[true_index]);
        LOGGER << "--------------------------------" << endl;
    }

    bool loop_break_indicator = false;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " "
            << excluded_SNP.size() << " " << backward_removed_SNP.size() << endl;

        for (int n : current_calculation_list)
            cohorts[n].calc_conditional_effects();

        inverse_var_meta(true);
        
        while (true) {
            // select maximal SNP
            double max_Zabs = sumstat_merge.col(2).maxCoeff(&max_screened_index);
            
            if (max_Zabs < 0 || sumstat_merge(max_screened_index, 3) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }
            
            true_index = screened_SNP[max_screened_index];
            string max_SNP_name = final_commonSNP[true_index];
            
            // calculate joint effects
            bool fail_flag = false;

            for (int n : current_calculation_list)
                append_row(cohorts[n].sumstat_candidate, cohorts[n].sumstat_screened.row(max_screened_index));
            
            for (int n : current_calculation_list) {
                auto &c = cohorts[n];
                c.calc_inner_product_with_SNP_list(candidate_SNP, true_index, window_size);
                c.calc_R_inv(if_fast_inv);

                if (c.calc_joint_effects(c.sumstat_candidate, true, iter_colinear_threshold)) {
                    // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                    fail_flag = true;
                    break;
                }

                if (c.R2 < (1+R2_incremental_threshold) * c.previous_R2) {
                    // LOGGER.w(1, "R2 increment unsatisfactory", max_SNP_name);
                    fail_flag = true;
                    break;
                }
            }

            if (fail_flag) {
                for (int n : current_calculation_list)
                    remove_row(cohorts[n].sumstat_candidate);

                sumstat_merge(max_screened_index, 2) = -1;
                continue;
            }

            // temporarily include new candidate SNP
            accept_SNP_as_candidate(true_index);
                
            // check joint model threshold
            inverse_var_meta_joint();

            if (sumstat_new_model_joint.col(3).bottomRows(candidate_SNP.size()-fixed_candidate_SNP_num).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed", max_SNP_name);
                    
                // LOGGER << "Added R vector cohort 1: " << c1.r.row(max_SNP_index) << endl;
                remove_new_colinear_SNP(true_index);

                int M = candidate_SNP.size();
                for (int n : current_calculation_list)
                    LOGGER << "Added diagonal value cohort " << n+1 << ": " << cohorts[n].R_inv_post(M-1, M-1) << endl;

                LOGGER << "Joint b: " << sumstat_new_model_joint(M-1, 0) << endl;
                LOGGER << "Joint se: " << sqrt(sumstat_new_model_joint(M-1, 1)) << endl;
                LOGGER << "Joint p-value: " << scientific << sumstat_new_model_joint(M-1, 3) << endl;
                
                for (int n : current_calculation_list)
                    LOGGER << "Adjusted R2 for cohort " << n+1 << ": " << fixed << cohorts[n].R2 << endl;

                break; 
            }

            // backward selection
            LOGGER << "Backward selection" << endl;

            for (int n : current_calculation_list)
                cohorts[n].sumstat_backward_new_model = cohorts[n].sumstat_candidate;

            vector<int> backward_new_model_candidate(candidate_SNP);

            for (int i = candidate_SNP.size()-1; i >= fixed_candidate_SNP_num; i--) {
                if (sumstat_new_model_joint(i, 3) > threshold) {
                    backward_new_model_candidate.erase(backward_new_model_candidate.begin() + i);
                    for (int n : current_calculation_list)
                        remove_row(cohorts[n].sumstat_backward_new_model, i);
                }
            }

            bool R2_backward_fail_flag = false;

            for (int n : current_calculation_list) {
                auto &c = cohorts[n];
                c.calc_R_inv_from_SNP_list(backward_new_model_candidate, window_size);
                c.calc_joint_effects(c.sumstat_backward_new_model, false);

                LOGGER << "Adjusted R2 after SNP elimination for cohort " << n+1 << ": " << fixed << c.R2 << endl; 
                if (c.R2 < (1+R2_incremental_threshold_backwards) * c.previous_R2) {
                    LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", max_SNP_name);
                    R2_backward_fail_flag = true;
                    break;
                }
            }

            if (R2_backward_fail_flag) {
                candidate_SNP.pop_back();

                for (int n : current_calculation_list) {
                    remove_row(cohorts[n].sumstat_candidate);
                    remove_column(cohorts[n].r);
                }

                sumstat_merge(max_screened_index, 2) = -1;
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", max_SNP_name);
            adjust_SNP_according_to_backward_selection();
            remove_new_colinear_SNP(true_index);
            break;
        }

        if (!loop_break_indicator)
            save_temp_model();

        LOGGER << "iter " << ++iter_num << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }
}
