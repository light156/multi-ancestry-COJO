#include "macojo.h"


void MACOJO::entry_function()
{   
    LOGGER.i(0, "Calculation started!");

    for (size_t n = 0; n < cohorts.size(); n++)
        current_list.push_back(n);

    if (params.if_joint_mode) {
        // --cojo-joint only works for single cohort, which is enforced in argument parsing
        candidate_SNP = screened_SNP; // screened SNP cannot be empty here, checked during reading files
        output_cojo_joint(params.output_name+".jma.cojo");
        LOGGER.i(0, "Calculation finished for joint analysis, program exit!");
        return;
    }  

    if (params.if_cond_mode) {
        // --cojo-cond only works for single cohort, which is enforced in argument parsing
        vector<string> temp;
        skim_file(params.cond_file, temp, false, true, false);
        
        for (const auto& SNP_name : temp) {
            auto iter = shared.goodSNP_index_map.find(SNP_name);
            if (iter == shared.goodSNP_index_map.end()){
                LOGGER.w(0, "conditional SNP not found in +[ " + params.cond_file + " ] and will be ignored", SNP_name);
                continue;
            }
            candidate_SNP.push_back(iter->second);
            screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), iter->second));
        }

        if (candidate_SNP.size() == 0)
            LOGGER.e(0, "None of the provided conditional SNPs exist in data");
    
        LOGGER.i(0, "effective conditional SNPs provided by the user", to_string(candidate_SNP.size()));

        // since only conditional analysis is performed, overriding main logic mode with effect size mode simplifies program flow
        params.if_gcta_COJO = (params.effect_size_mode == "GCTA");
        params.if_remove_NA = (params.effect_size_mode == "removeNA");

        for (int index : candidate_SNP)
            cohorts[0].update_r(screened_SNP, index);

        if (!cohorts[0].calc_R_inv_from_SNP_list(candidate_SNP, params.if_gcta_COJO, params.if_remove_NA))
            LOGGER.w(0, "Colinearity exceeds threshold when calculating R inverse for conditional SNPs, make sure this is what you want"); 
        
        cohorts[0].save_temp_model();

        output_cojo_cond(params.output_name+".cma.cojo");
        LOGGER.i(0, "Calculation finished for conditional analysis, program exit!");
        return;
    }

    // COJO stepwise iterative selection, first deal with fixed candidate SNPs if provided
    if (!params.fixedSNP_file.empty()) {
        vector<string> temp;
        skim_file(params.fixedSNP_file, temp, false, true, false);

        for (const auto& SNP_name : temp) {
            auto iter = shared.goodSNP_index_map.find(SNP_name);
            if (iter == shared.goodSNP_index_map.end()){
                LOGGER.w(0, "fixed candidate SNP not found in +[ " + params.fixedSNP_file + " ] and will be ignored", SNP_name);
                continue;
            }
            candidate_SNP.push_back(iter->second);
            screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), iter->second));
        }

        if (candidate_SNP.size() == 0)
            LOGGER.e(0, "None of the provided fixed candidate SNPs exist in data");

        LOGGER.i(0, "effective fixed candidate SNPs provided by the user", to_string(candidate_SNP.size()));

        for (size_t n = 0; n < cohorts.size(); n++) {
            auto &c = cohorts[n];

            for (int index : candidate_SNP)
                c.update_r(screened_SNP, index);

            if (!c.calc_R_inv_from_SNP_list(candidate_SNP, params.if_gcta_COJO, params.if_remove_NA)) {
                LOGGER << "R matrix for Cohort " << n+1 << ":" << endl << c.R_post << endl;
                LOGGER.e(0, "Colinearity exceeds threshold when calculating R inverse for fixed candidate SNPs, please check"); 
            }

            if (!c.calc_joint_effects(params.if_gcta_COJO)) {
                LOGGER << "R matrix for Cohort " << n+1 << ":" << endl << c.R_post << endl;
                LOGGER.e(0, "Joint se too small for fixed candidate SNPs, please check");
            }

            c.save_temp_model();
        }
    }

    fixed_candidate_SNP_num = candidate_SNP.size();
    main_loop();
    output_cojo_joint(params.output_name + ".jma.cojo");

    if (cohorts.size() > 1 && params.if_MDISA) {
        current_list.resize(1);

        // backup SNPs for MDISA, Manc-COJO candidate SNPs will be fixed
        fixed_candidate_SNP_num = candidate_SNP.size();
        candidate_SNP_backup = candidate_SNP;

        for (size_t temp_index = 0; temp_index < cohorts.size(); temp_index++) {
            LOGGER << "MDISA for Cohort " << temp_index+1 << endl;
            current_list[0] = temp_index;
            
            candidate_SNP = candidate_SNP_backup;
            screened_SNP.insert(screened_SNP.end(), removed_SNP.begin(), removed_SNP.end());
            vector<int>().swap(removed_SNP);

            main_loop();
            output_cojo_joint(params.output_name + ".MDISA.cohort" + to_string(temp_index + 1) + ".jma.cojo");
        }
    }

    LOGGER.i(0, "All calculations finished!");
}


void MACOJO::main_loop() 
{   
    int min_pC_index, max_pJ_index;
    double min_pC, max_pJ;

    string current_SNP_name;
    bool loop_break_indicator = false;
    bool skip_flag, backward_success_flag;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num < params.max_iter_num && screened_SNP.size() > 0) {
        auto start = chrono::steady_clock::now();

        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " " << removed_SNP.size() << endl;

        inverse_var_meta_conditional();

        if (candidate_SNP.size() == 0) {
            LOGGER.i(0, "No candidate SNPs, using the most significant SNP as the first candidate");

            min_pC = erfc(abs_zC.maxCoeff(&min_pC_index) / sqrt(2));
            if (min_pC > params.p)  {// erfc(-1) = 1.8427 > 1 so don't need to worry about abs_zC being -1
                LOGGER.w(0, "Calculation finished, Input data has no significant SNPs");
                break;
            };

            current_SNP_name = shared.SNP_ref[min_pC_index];
            LOGGER.i(0, "First SNP", current_SNP_name);
            LOGGER << "p-value: " << scientific << min_pC << fixed << endl;
            LOGGER << "--------------------------------" << endl;

            candidate_SNP.push_back(min_pC_index);
            screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));

            for (int n : current_list) {
                cohorts[n].sumstat_candidate = cohorts[n].sumstat.row(min_pC_index);
                cohorts[n].update_r(screened_SNP, min_pC_index);
                cohorts[n].save_temp_model();
            }

            iter_num++;
            continue;
        }
        
        while (true) {
            skip_flag = false;
            backward_success_flag = true;

            // select minimal conditional SNP
            min_pC = erfc(abs_zC.maxCoeff(&min_pC_index) / sqrt(2));
            current_SNP_name = shared.SNP_ref[min_pC_index];

            if (min_pC > params.p) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                LOGGER << current_SNP_name << " " << scientific << min_pC << fixed << endl;
                loop_break_indicator = true;
                break;
            }
            
            // calculate R inverse for the new candidate SNP
            for (int n : current_list) {
                if (!cohorts[n].calc_R_inv_forward(min_pC_index)) {
                    // LOGGER.w(1, "skipped, colinearity exceeds threshold in Cohort " + to_string(n+1), current_SNP_name);
                    skip_flag = true;
                    break;
                }
            }
            
            if (skip_flag) {
                abs_zC(min_pC_index) = -1;
                continue;
            }
    
            // temporarily update sumstat_candidate for calculating joint effects
            for (int n : current_list)
                append_row(cohorts[n].sumstat_candidate, cohorts[n].sumstat.row(min_pC_index));

            for (int n : current_list) {  
                if (!cohorts[n].calc_joint_effects(params.if_gcta_COJO)) {
                    // LOGGER.w(1, "skipped, joint se too small in Cohort " + to_string(n+1), current_SNP_name);
                    skip_flag = true;
                    break;
                }
            }

            // remove this SNP if joint se too small
            if (skip_flag) {
                for (int n : current_list)
                    remove_row(cohorts[n].sumstat_candidate);

                removed_SNP.push_back(min_pC_index);
                screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));

                abs_zC(min_pC_index) = -1;
                continue;
            }

            // check R2 increment, which does not need to be done in gcta-COJO
            if (!params.if_gcta_COJO)  {
                for (int n : current_list) {
                    if (cohorts[n].R2 < (1+params.R2_threshold) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "skipped, R2 increment lower than threshold", current_SNP_name);
                        skip_flag = true;
                    }
                }
                
                if (skip_flag) {
                    for (int n : current_list)
                        remove_row(cohorts[n].sumstat_candidate);

                    abs_zC(min_pC_index) = -1;
                    continue;
                }
            }
            
            LOGGER << "Conditional min pC: " << scientific << min_pC << fixed << endl;

            // calculate joint p-value and decide whether to accept this SNP
            inverse_var_meta_joint();

            max_pJ = erfc(abs_zJ.minCoeff(&max_pJ_index) / sqrt(2));
            LOGGER << "Joint b, se, max pJ: " 
                    << bJ(max_pJ_index) << " " << sqrt(se2J(max_pJ_index)) << " " << scientific << max_pJ << fixed << endl;

            // directly accept this SNP and next iteration
            if (max_pJ <= params.p) {
                LOGGER.i(0, "New candidate SNP accepted", current_SNP_name);

                candidate_SNP.push_back(min_pC_index);
                screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));

                for (int n : current_list) {
                    cohorts[n].update_r(screened_SNP, min_pC_index);
                    cohorts[n].save_temp_model();
                }

                break;
            } 
            
            // trivial case for removing current SNP
            if (max_pJ_index < fixed_candidate_SNP_num || max_pJ_index == abs_zJ.rows()-1) {
                LOGGER.i(0, "removed, because pJ exceeds p-value threshold", current_SNP_name);
                for (int n : current_list)
                    remove_row(cohorts[n].sumstat_candidate);

                removed_SNP.push_back(min_pC_index);
                screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));

                abs_zC(min_pC_index) = -1;
                continue;
            } 

            // backward selection must happen now, save current model as backup, then accept the new candidate SNP
            candidate_SNP_backup = candidate_SNP;
            removed_SNP_backup = removed_SNP;

            for (int n : current_list)
                cohorts[n].save_state();

            LOGGER.i(0, "New candidate SNP added, start backward selection", current_SNP_name);

            candidate_SNP.push_back(min_pC_index);
            screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));

            for (int n : current_list) {
                cohorts[n].update_r(screened_SNP, min_pC_index);
                cohorts[n].save_temp_model();
            }

            // backward selection
            while (true) {
                LOGGER.i(0, "Candidate SNP removed", shared.SNP_ref[candidate_SNP[max_pJ_index]]);
                
                // remove r_matrix column and sumstat_candidate row
                for (int n : current_list) {
                    remove_row(cohorts[n].sumstat_candidate, max_pJ_index);
                    remove_column(cohorts[n].r, max_pJ_index);

                    if (params.if_gcta_COJO)
                        remove_column(cohorts[n].r_gcta, max_pJ_index);
                }

                removed_SNP.push_back(candidate_SNP[max_pJ_index]);
                candidate_SNP.erase(candidate_SNP.begin() + max_pJ_index);

                if (candidate_SNP.size() == 0) {
                    LOGGER.i(0, "Backward selection successful, no candidate SNPs left");
                    break;
                }

                for (int n : current_list) {
                    // do not need to check colinearity because it is guaranteed in forward step
                    cohorts[n].calc_R_inv_backward(max_pJ_index);

                    if (!cohorts[n].calc_joint_effects(params.if_gcta_COJO)) {
                        // don't think this can happen, but just in case
                        LOGGER.i(0, "Backward selection failed, joint se too small after removing " + current_SNP_name);
                        backward_success_flag = false;
                        break;
                    }
                }

                if (!backward_success_flag) break;

                inverse_var_meta_joint();

                max_pJ = erfc(abs_zJ.minCoeff(&max_pJ_index) / sqrt(2));
                LOGGER << "Joint b, se, max pJ: " 
                        << bJ(max_pJ_index) << " " << sqrt(se2J(max_pJ_index)) << " " << scientific << max_pJ << fixed << endl;

                for (int n : current_list)
                    cohorts[n].save_temp_model();

                if (max_pJ <= params.p) break;

                if (max_pJ_index < fixed_candidate_SNP_num) {
                    LOGGER.i(0, "Backward selection failed, fixed candidate SNP exceeds p-value threshold");
                    backward_success_flag = false;
                    break;
                }
            } 

            // no need to save temp model here because the model will start from zero
            if (candidate_SNP.size() == 0) break;
            
            // check R2 increment after backward selection, which does not need to be done in gcta-COJO
            if (backward_success_flag && !params.if_gcta_COJO) {
                for (int n : current_list) {
                    if (cohorts[n].R2 < (1+params.R2back_threshold) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "Backward selection failed, adjusted R2 lower than backward threshold", current_SNP_name);
                        backward_success_flag = false;
                    }
                }
            }

            // restore the last successful model before the newest candidate SNP is added
            // the last row of sumstat_candidate in the restored state needs to be removed
            if (!backward_success_flag) {
                LOGGER.i(0, "Restore to the last successful model before adding " + current_SNP_name);
                candidate_SNP = candidate_SNP_backup;
                removed_SNP = removed_SNP_backup;
                for (int n : current_list) {
                    cohorts[n].restore_state();
                    remove_row(cohorts[n].sumstat_candidate); 
                }

                abs_zC(min_pC_index) = -1;
                continue;
            }

            LOGGER.i(0, "Backward selection successful!");
            break;
        }

        auto end = chrono::steady_clock::now();
        LOGGER << "iter " << iter_num++ << " finished, takes " 
                << chrono::duration<double>(end-start).count() << " seconds" << endl;
        LOGGER << "--------------------------------" << endl;
    }
}
