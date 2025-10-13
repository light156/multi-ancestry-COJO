#include "macojo.h"


void MACOJO::entry_function(string savename)
{   
    LOGGER.i(0, "Calculation started!");

    for (int n = 0; n < cohorts.size(); n++)
        current_calculation_list.push_back(n);

    initialize_main_loop();

    if (if_cojo_joint){
        LOGGER.i(0, "COJO-Joint analysis only, no iterative selection");

        if (fixed_candidate_SNP_num == 0)
            LOGGER.w(0, "No candidate SNPs provided for COJO-Joint analysis, program exit!");
        else {
            output_results_to_file(savename+".jma.cojo.fixedSNP");
            LOGGER.i(0, "Calculation finished with fixed candidate SNPs, program exit!");
        }
        return;
    }

    main_loop();
    output_results_to_file(savename+".jma.cojo");

    if (cohorts.size() > 1 && !if_skip_MDISA) {
        // backup SNPs for MDISA
        fixed_candidate_SNP_num = candidate_SNP.size();
        candidate_SNP_backup = candidate_SNP;
        screened_SNP_backup = screened_SNP;
        removed_SNP_backup = removed_SNP;

        current_calculation_list.resize(1);

        for (int temp_index = 0; temp_index < cohorts.size(); temp_index++) {
            LOGGER << "MDISA for Cohort " << temp_index+1 << endl;
            current_calculation_list[0] = temp_index;
            
            candidate_SNP = candidate_SNP_backup;
            screened_SNP = screened_SNP_backup;
            removed_SNP = removed_SNP_backup;

            initialize_MDISA();
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

    if (fixed_candidate_SNP.size() == 0)
        return;
    
    list<int> candidate_SNP_temp;

    for (const auto &SNP_name : fixed_candidate_SNP) {
        int single_index = distance(final_commonSNP.begin(), find(final_commonSNP.begin(), final_commonSNP.end(), SNP_name));
        if (single_index == final_commonSNP.size())
            LOGGER.w(0, "Fixed candidate SNP not found in data and will be ignored", SNP_name);
        else
            candidate_SNP_temp.push_back(single_index);
    }

    if (candidate_SNP_temp.size() == 0) {
        LOGGER.w(0, "None of the fixed candidate SNPs provided by the user exist in data, program continues without fixed candidate SNPs");
        return;
    }

    candidate_SNP_temp.unique(); // currently candidate_SNP_temp may have duplicates
    candidate_SNP_temp.sort(greater<int>()); // sort in decreasing order to avoid index change when erasing from screened_SNP

    for (int single_index : candidate_SNP_temp) {
        if (!if_cojo_joint) {
            accept_SNP_as_candidate(single_index);
        } else {
            candidate_SNP.push_back(screened_SNP[single_index]);
            for (int n : current_calculation_list) 
                append_row(cohorts[n].sumstat_candidate, cohorts[n].sumstat_screened.row(single_index));
        }
    }

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];

        if (!c.calc_R_inv_from_SNP_list(candidate_SNP, c.sumstat_candidate))
            LOGGER.e(0, "Colinearity exceeds threshold when calculating R inverse for fixed candidate SNPs, please check");
        
        c.calc_joint_effects();

        if (c.beta_var.minCoeff() < 1e-30)
            LOGGER.e(0, "Joint se too small for fixed candidate SNPs, please check");

        c.save_temp_model();
    }
    
    LOGGER.i(0, "effective fixed candidate SNPs provided by the user", to_string(fixed_candidate_SNP_num));
}


void MACOJO::main_loop() 
{   
    int max_screened_index, max_joint_index;
    string max_SNP_name;

    bool loop_break_indicator = false;
    bool skip_flag, backward_success_flag = false;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " " << removed_SNP.size() << endl;

        if (candidate_SNP.size() == 0) {
            LOGGER.i(0, "No candidate SNPs, using the most significant SNP as the first candidate");

            inverse_var_meta_init();

            sumstat_merge.col(2).maxCoeff(&max_screened_index);
            if (sumstat_merge(max_screened_index, 3) > threshold)
                LOGGER.e(0, "Input data has no significant SNPs");
            
            max_SNP_name = final_commonSNP[screened_SNP[max_screened_index]];
            
            accept_SNP_as_candidate(max_screened_index);
            LOGGER.i(0, "First SNP", max_SNP_name);
            LOGGER << "p-value: " << scientific << sumstat_merge(max_screened_index, 3) << fixed << endl;
            LOGGER << "--------------------------------" << endl;

            for (int n : current_calculation_list)
                cohorts[n].save_temp_model();

            iter_num++;
            continue;
        }

        for (int n : current_calculation_list)
            cohorts[n].calc_conditional_effects();

        inverse_var_meta_conditional();
        
        while (true) {
            skip_flag = false;
            backward_success_flag = false;

            // select maximal SNP
            double max_Zabs = sumstat_merge.col(2).maxCoeff(&max_screened_index);
            max_SNP_name = final_commonSNP[screened_SNP[max_screened_index]];

            if (max_Zabs < 0 || sumstat_merge(max_screened_index, 3) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                LOGGER << max_SNP_name << " " << scientific << sumstat_merge(max_screened_index, 3) << fixed << endl;
                loop_break_indicator = true;
                break;
            }

            for (int n : current_calculation_list) {
                auto &c = cohorts[n];
                
                if (c.sumstat_screened(max_screened_index, 6) > 1e10) {
                    // LOGGER.w(1, "skipped, adjusted N too large", max_SNP_name);
                    skip_flag = true;
                    break;
                }

                c.r_temp_vec = c.r.row(max_screened_index).transpose();

                // gcta-COJO colinearity check
                if (if_gcta_COJO && (c.r_temp_vec.transpose() * c.R_inv_pre * c.r_temp_vec > colinear_threshold)) {
                    // LOGGER.w(1, "skipped, colinear with existing candidate SNPs", max_SNP_name);
                    skip_flag = true;
                    break;
                }

                // our own colinearity check
                if (!if_gcta_COJO && (c.r_temp_vec.cwiseAbs2().maxCoeff() > colinear_threshold)) {
                    // LOGGER.w(1, "skipped, colinear with existing candidate SNPs", max_SNP_name);
                    skip_flag = true;
                    break;
                }
                
                if (!c.calc_R_inv()) {
                    // LOGGER.w(1, "skipped, NA produced during R inverse calculation", max_SNP_name);
                    skip_flag = true;
                    break;
                }
                
                if (if_gcta_COJO) 
                    c.calc_R_inv_gcta(max_screened_index);
                
                c.sumstat_new_model = c.sumstat_candidate;
                append_row(c.sumstat_new_model, c.sumstat_screened.row(max_screened_index));
                if (!c.calc_joint_effects()) {
                    LOGGER.w(1, "skipped, joint se too small", max_SNP_name);
                    skip_flag = true;
                    break;
                }
            }

            // check R2 increment after forward selection, which does not need to be done in gcta-COJO
            if (!skip_flag && !if_gcta_COJO) {
                for (int n : current_calculation_list) {
                    if (cohorts[n].R2 < (1+R2_incremental_threshold) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "skipped, R2 increment lower than threshold", max_SNP_name);
                        skip_flag = true;
                    }
                }
            }
                        
            if (skip_flag) {
                sumstat_merge(max_screened_index, 2) = -1;
                continue;
            }
    
            LOGGER << "Conditional: " << sumstat_merge(max_screened_index, 0) << " "
                    << sqrt(sumstat_merge(max_screened_index, 1)) << " "
                    << scientific << sumstat_merge(max_screened_index, 3) << fixed << endl;

            inverse_var_meta_joint();
            double max_joint_pJ = sumstat_merge_new_model.col(3).maxCoeff(&max_joint_index);
            LOGGER << "Joint max pJ: " << sumstat_merge_new_model(max_joint_index, 0) << " "
                    << sqrt(sumstat_merge_new_model(max_joint_index, 1)) << " "
                    << scientific << max_joint_pJ << fixed << endl;
        
            if (max_joint_pJ <= threshold) {
                LOGGER.i(0, "New candidate SNP accepted", max_SNP_name);
                accept_SNP_as_candidate(max_screened_index);
                for (int n : current_calculation_list)
                    cohorts[n].save_temp_model();
                break;
            }

            backward_success_flag = true;
            candidate_SNP_backup = candidate_SNP;
            removed_SNP_backup = removed_SNP;
            for (int n : current_calculation_list)
                cohorts[n].save_state();

            LOGGER.i(0, "New candidate SNP added, start backward selection", max_SNP_name);
            accept_SNP_as_candidate(max_screened_index);
            for (int n : current_calculation_list)
                cohorts[n].save_temp_model();

            // backward selection
            while (max_joint_pJ > threshold) {
                if (max_joint_index < fixed_candidate_SNP_num) {
                    LOGGER.i(0, "Backward selection failed, because fixed candidate SNP has high joint p-value");
                    backward_success_flag = false;
                    break;
                }

                LOGGER.i(0, "Removed candidate SNP", final_commonSNP[candidate_SNP[max_joint_index]]);
                removed_SNP.push_back(candidate_SNP[max_joint_index]);
                candidate_SNP.erase(candidate_SNP.begin() + max_joint_index);

                for (int n : current_calculation_list) {
                    auto &c = cohorts[n];
                    append_row(c.sumstat_removed, c.sumstat_candidate.row(max_joint_index));
                    remove_row(c.sumstat_candidate, max_joint_index);
                    remove_column(c.r, max_joint_index);
                    if (if_gcta_COJO)
                        remove_column(c.r_gcta, max_joint_index);
                }

                if (candidate_SNP.size() == 0) break;

                for (int n : current_calculation_list) {
                    auto &c = cohorts[n];
                    c.calc_R_inv_backward(max_joint_index); // do not need to check colinearity

                    remove_row(c.sumstat_new_model, max_joint_index);
                    if (!c.calc_joint_effects()) {
                        LOGGER.i(0, "Backward selection failed, joint se too small after removing " + max_SNP_name);
                        backward_success_flag = false;
                        break;
                    }
                }

                if (!backward_success_flag) break;

                inverse_var_meta_joint();
                max_joint_pJ = sumstat_merge_new_model.col(3).maxCoeff(&max_joint_index);
                LOGGER << "Joint max pJ: " << sumstat_merge_new_model(max_joint_index, 0) << " "
                        << sqrt(sumstat_merge_new_model(max_joint_index, 1)) << " "
                        << scientific << max_joint_pJ << fixed << endl;

                for (int n : current_calculation_list)
                    cohorts[n].save_temp_model();
            }

            if (candidate_SNP.size() == 0) {
                LOGGER.i(0, "Backward selection successful, no candidate SNPs left");
                break;
            }
            
            if (backward_success_flag && !if_gcta_COJO) {
                for (int n : current_calculation_list) {
                    if (cohorts[n].R2 < (1+R2_incremental_threshold_backwards) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "Backward selection failed, adjusted R2 lower than backward threshold", max_SNP_name);
                        backward_success_flag = false;
                    }
                }
            }

            // restore the last successful model before the newest candidate SNP is added
            if (!backward_success_flag) {
                LOGGER.i(0, "Restore to the last successful model before adding " + max_SNP_name);
                candidate_SNP = candidate_SNP_backup;
                removed_SNP = removed_SNP_backup;
                for (int n : current_calculation_list)
                    cohorts[n].restore_state();
            } else
                LOGGER.i(0, "Backward selection successful");
            
            break;
        }

        LOGGER << "iter " << iter_num++ << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }
}
