#include "macojo.h"


void MACOJO::entry_function(string savename)
{   
    LOGGER.i(0, "Calculation started!");

    for (int n = 0; n < cohorts.size(); n++)
        current_calculation_list.push_back(n);

    initialize_main_loop();

    if (params.if_cojo_joint){
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

    if (cohorts.size() > 1 && !params.if_skip_MDISA) {
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
    screened_SNP.resize(shared.commonSNP_total_num);
    for (int i = 0; i < shared.commonSNP_total_num; i++)
        screened_SNP[i] = i;

    if (fixed_candidate_SNP.size() == 0)
        return;
    
    LOGGER.i(0, "The user has provided fixed candidate SNPs for analysis");

    list<int> candidate_SNP_temp;
    for (const auto &SNP_name : fixed_candidate_SNP) {
        int single_index = distance(shared.final_commonSNP.begin(), find(shared.final_commonSNP.begin(), shared.final_commonSNP.end(), SNP_name));
        if (single_index == shared.final_commonSNP.size())
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
        for (int n : current_calculation_list) 
            append_row(cohorts[n].sumstat_candidate, cohorts[n].sumstat_screened.row(single_index));

        if (params.if_cojo_joint)
            candidate_SNP.push_back(screened_SNP[single_index]);
        else
            accept_SNP_as_candidate(single_index);
    }

    for (int n : current_calculation_list) {
        auto &c = cohorts[n];

        if (!c.calc_R_inv_from_SNP_list(candidate_SNP, c.sumstat_candidate)) {
            LOGGER << "R matrix for Cohort " << n+1 << ":" << endl << c.R_post << endl;
            LOGGER.e(0, "Colinearity exceeds threshold when calculating R inverse for fixed candidate SNPs, please check");
        }

        if (!c.calc_joint_effects()) {
            LOGGER << "R matrix for Cohort " << n+1 << ":" << endl << c.R_post << endl;
            LOGGER.e(0, "Joint se too small for fixed candidate SNPs, please check");
        }

        c.save_temp_model();
    }
    
    LOGGER.i(0, "effective fixed candidate SNPs provided by the user", to_string(fixed_candidate_SNP_num));
}


void MACOJO::main_loop() 
{   
    int min_screened_index, max_joint_index;
    double min_pC, max_pJ;

    string current_SNP_name;
    bool loop_break_indicator = false;
    bool skip_flag, backward_success_flag;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num < params.max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " " << removed_SNP.size() << endl;

        if (candidate_SNP.size() == 0) {
            LOGGER.i(0, "No candidate SNPs, using the most significant SNP as the first candidate");

            inverse_var_meta_init();
            min_pC = erfc(abs_zC.maxCoeff(&min_screened_index) / sqrt(2)) / 2;
            if (min_pC > params.threshold)
                LOGGER.e(0, "Input data has no significant SNPs");

            current_SNP_name = shared.final_commonSNP[screened_SNP[min_screened_index]];
            LOGGER.i(0, "First SNP", current_SNP_name);
            LOGGER << "p-value: " << scientific << min_pC << fixed << endl;
            LOGGER << "--------------------------------" << endl;

            for (int n : current_calculation_list) 
                cohorts[n].sumstat_candidate = cohorts[n].sumstat_screened.row(min_screened_index);

            accept_SNP_as_candidate(min_screened_index);
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
            backward_success_flag = true;

            // select minimal conditional SNP
            min_pC = erfc(abs_zC.maxCoeff(&min_screened_index) / sqrt(2)) / 2;
            current_SNP_name = shared.final_commonSNP[screened_SNP[min_screened_index]];

            if (abs_zC(min_screened_index) < 0 || min_pC > params.threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                LOGGER << current_SNP_name << " " << scientific << min_pC << fixed << endl;
                loop_break_indicator = true;
                break;
            }

            LOGGER.i(0, "Screened SNP with the most significant conditional pC", current_SNP_name);

            // temporarily update sumstat_candidate for calculating joint effects
            for (int n : current_calculation_list) 
                append_row(cohorts[n].sumstat_candidate, cohorts[n].sumstat_screened.row(min_screened_index));

            // calculate joint effects and check colinearity for each cohort
            for (int n : current_calculation_list) {
                auto &c = cohorts[n];
                
                if (c.sumstat_screened(min_screened_index, 6) > 1e10) {
                    // LOGGER.w(1, "skipped, adjusted N too large", current_SNP_name);
                    skip_flag = true;
                    break;
                }

                c.r_temp_vec = c.r.row(min_screened_index).transpose();
                
                // colinearity check
                // c.r_temp_vec.cwiseAbs2().maxCoeff() > params.colinear_threshold
                if (c.r_temp_vec.transpose() * c.R_inv_pre * c.r_temp_vec > params.colinear_threshold) {
                    // LOGGER.w(1, "skipped, colinear with existing candidate SNPs", current_SNP_name);
                    skip_flag = true;
                    break;
                }
                
                if (!c.calc_R_inv()) {
                    // LOGGER.w(1, "skipped, NA produced during R inverse calculation", current_SNP_name);
                    skip_flag = true;
                    break;
                }
                
                if (params.if_gcta_COJO) 
                    c.calc_R_inv_gcta(min_screened_index);
                
                if (!c.calc_joint_effects()) {
                    LOGGER.w(1, "skipped, joint se too small", current_SNP_name);
                    skip_flag = true;
                    break;
                }
            }

            // check R2 increment after forward selection, which does not need to be done in gcta-COJO
            if (!skip_flag && !params.if_gcta_COJO) {
                for (int n : current_calculation_list) {
                    if (cohorts[n].R2 < (1+params.R2_incremental_threshold) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "skipped, R2 increment lower than threshold", current_SNP_name);
                        skip_flag = true;
                    }
                }
            }
            
            // if anything wrong, revert sumstat_candidate and try the next SNP
            if (skip_flag) {
                for (int n : current_calculation_list)
                    remove_row(cohorts[n].sumstat_candidate);
                
                abs_zC(min_screened_index) = -1;
                continue;
            }
    
            LOGGER << "Conditional min pC: " << scientific << min_pC << fixed << endl;

            // calculate joint p-value and decide whether to accept this SNP
            inverse_var_meta_joint();

            max_pJ = erfc(abs_zJ.minCoeff(&max_joint_index) / sqrt(2)) / 2;
            LOGGER << "Joint b, se, max pJ: " 
                    << bJ(max_joint_index) << " " << sqrt(se2J(max_joint_index)) << " " << scientific << max_pJ << fixed << endl;
            
            // directly accept this SNP and next iteration
            if (max_pJ <= params.threshold) {
                LOGGER.i(0, "New candidate SNP accepted", current_SNP_name);
                accept_SNP_as_candidate(min_screened_index);
                for (int n : current_calculation_list)
                    cohorts[n].save_temp_model();

                break;
            } 
            
            // trivial case for removing this SNP
            if (max_joint_index < fixed_candidate_SNP_num) {
                LOGGER.i(0, "skipped, fixed candidate SNP exceeds p-value threshold", current_SNP_name);
                for (int n : current_calculation_list)
                    remove_row(cohorts[n].sumstat_candidate);

                remove_row(abs_zC, min_screened_index);
                remove_SNP_from_screened(min_screened_index);
                continue;
            } 

            // trivial case for removing this SNP
            if (max_joint_index == abs_zJ.rows()-1) {
                LOGGER.i(0, "skipped, it has the highest joint p-value", current_SNP_name);
                for (int n : current_calculation_list)
                    remove_row(cohorts[n].sumstat_candidate);

                remove_row(abs_zC, min_screened_index);
                remove_SNP_from_screened(min_screened_index);
                continue;
            }

            // backward selection must happen now, save current model as backup
            candidate_SNP_backup = candidate_SNP;
            screened_SNP_backup = screened_SNP;
            removed_SNP_backup = removed_SNP;
            for (int n : current_calculation_list)
                cohorts[n].save_state();

            LOGGER.i(0, "New candidate SNP added, start backward selection", current_SNP_name);

            accept_SNP_as_candidate(min_screened_index);
            for (int n : current_calculation_list)
                cohorts[n].save_temp_model();

            // backward selection
            while (true) {
                LOGGER.i(0, "Candidate SNP removed", shared.final_commonSNP[candidate_SNP[max_joint_index]]);
                remove_SNP_from_candidate(max_joint_index);

                if (candidate_SNP.size() == 0) break;

                for (int n : current_calculation_list) {
                    auto &c = cohorts[n];
                    // do not need to check colinearity because it is guaranteed in forward step
                    calc_R_inverse_backward(c.R_pre, c.R_inv_pre, max_joint_index, c.R_post, c.R_inv_post, params.if_fast_inv);
                    if (params.if_gcta_COJO)
                        calc_R_inverse_backward(c.R_pre_gcta, c.R_inv_pre_gcta, max_joint_index, 
                                                c.R_post_gcta, c.R_inv_post_gcta, params.if_fast_inv);

                    if (!c.calc_joint_effects()) {
                        LOGGER.i(0, "Backward selection failed, joint se too small after removing " + current_SNP_name);
                        backward_success_flag = false;
                        break;
                    }
                }

                if (!backward_success_flag) break;

                inverse_var_meta_joint();

                max_pJ = erfc(abs_zJ.minCoeff(&max_joint_index) / sqrt(2)) / 2;
                LOGGER << "Joint b, se, max pJ: " 
                        << bJ(max_joint_index) << " " << sqrt(se2J(max_joint_index)) << " " << scientific << max_pJ << fixed << endl;

                for (int n : current_calculation_list)
                    cohorts[n].save_temp_model();

                if (max_pJ <= params.threshold) break;

                if (max_joint_index < fixed_candidate_SNP_num) {
                    LOGGER.i(0, "Backward selection failed, fixed candidate SNP exceeds p-value threshold");
                    backward_success_flag = false;
                    break;
                }
            } 

            // no need to save temp model here because the model will start from zero
            if (candidate_SNP.size() == 0) {
                LOGGER.i(0, "Backward selection successful, no candidate SNPs left");
                break;
            }
            
            // check R2 increment after backward selection, which does not need to be done in gcta-COJO
            if (backward_success_flag && !params.if_gcta_COJO) {
                for (int n : current_calculation_list) {
                    if (cohorts[n].R2 < (1+params.R2_incremental_threshold_backwards) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "Backward selection failed, adjusted R2 lower than backward threshold", current_SNP_name);
                        backward_success_flag = false;
                    }
                }
            }

            // restore the last successful model before the newest candidate SNP is added
            if (!backward_success_flag) {
                LOGGER.i(0, "Restore to the last successful model before adding " + current_SNP_name);
                candidate_SNP = candidate_SNP_backup;
                screened_SNP = screened_SNP_backup;
                removed_SNP = removed_SNP_backup;
                for (int n : current_calculation_list) {
                    cohorts[n].restore_state();
                    remove_row(cohorts[n].sumstat_candidate);
                }
                
                abs_zC(min_screened_index) = -1;
                continue;
            }

            LOGGER.i(0, "Backward selection successful!");
            break;
        }

        LOGGER << "iter " << iter_num++ << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }
}
