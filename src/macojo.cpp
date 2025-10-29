#include "macojo.h"


void MACOJO::entry_function()
{   
    LOGGER.i(0, "Calculation started!");

    for (size_t n = 0; n < cohorts.size(); n++)
        current_list.push_back(n);

    initialize_from_fixed_candidate();
    if (params.if_cojo_joint) return;

    fixed_candidate_SNP_num = candidate_SNP.size();
    main_loop();
    output_results_to_file(output_name+".jma.cojo");

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
            output_results_to_file(output_name+".MDISA.cohort" + to_string(temp_index+1) + ".jma.cojo");
        }
    }

    LOGGER.i(0, "All calculations finished!");
}


void MACOJO::initialize_from_fixed_candidate() 
{   
    if (params.if_cojo_joint) {
        // --cojo-joint only works for single cohort, which is enforced in argument parsing
        candidate_SNP = screened_SNP; // screened SNP cannot be empty here, checked during reading files
        bool flag = true;

        if (!cohorts[0].calc_R_inv_from_SNP_list(candidate_SNP)) {
            LOGGER.i(0, "Colinearity exceeds threshold when calculating R inverse, make sure this is what you want");
            flag = false;
        }

        cohorts[0].sumstat_candidate = cohorts[0].sumstat(candidate_SNP, all);

        if (!cohorts[0].calc_joint_effects()) {
            LOGGER.i(0, "Joint se too small, make sure this is what you want");
            flag = false;
        }

        cohorts[0].save_temp_model();

        if (flag)
            output_results_to_file(output_name+".jma.cojo");
        else {   
            LOGGER << "R matrix:" << endl << cohorts[0].R_post << endl;
            output_results_to_file(output_name+".jma.cojo.warning");
        }
            
        LOGGER.i(0, "Calculation finished for COJO-Joint analysis, program exit!");
        return;
    }  

    if (fixed_candidate_SNP.size() == 0) return;
    
    list<int> candidate_SNP_temp;
    for (const auto &SNP_name : fixed_candidate_SNP) {
        auto iter = shared.goodSNP_index_map.find(SNP_name);
        if (iter == shared.goodSNP_index_map.end()){
            LOGGER.w(0, "Fixed candidate SNP not found in data and will be ignored", SNP_name);
            continue;
        }
        candidate_SNP_temp.push_back(iter->second);
    }

    if (candidate_SNP_temp.size() == 0)
        LOGGER.e(0, "None of the fixed candidate SNPs provided by the user exist in data");

    candidate_SNP_temp.unique();
    LOGGER.i(0, "effective fixed candidate SNPs provided by the user", to_string(candidate_SNP_temp.size()));

    for (int index : candidate_SNP_temp) {
        for (auto& c : cohorts) 
            append_row(c.sumstat_candidate, c.sumstat.row(index));
        
        accept_SNP_as_candidate(index);
    }
    
    for (size_t n = 0; n < cohorts.size(); n++) {
        auto &c = cohorts[n];

        if (!c.calc_R_inv_from_SNP_list(candidate_SNP)) {
            LOGGER << "R matrix for Cohort " << n+1 << ":" << endl << c.R_post << endl;
            LOGGER.e(0, "Colinearity exceeds threshold when calculating R inverse for fixed candidate SNPs, please check");
        }

        if (!c.calc_joint_effects()) {
            LOGGER << "R matrix for Cohort " << n+1 << ":" << endl << c.R_post << endl;
            LOGGER.e(0, "Joint se too small for fixed candidate SNPs, please check");
        }

        c.save_temp_model();
    }
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

        if (candidate_SNP.size() == 0) {
            LOGGER.i(0, "No candidate SNPs, using the most significant SNP as the first candidate");

            inverse_var_meta_conditional();
            min_pC = erfc(abs_zC.maxCoeff(&min_pC_index) / sqrt(2));
            if (min_pC > params.p)  // erfc(-1) = 1.8427 > 1 so don't need to worry about abs_zC being -1
                LOGGER.e(0, "Input data has no significant SNPs");

            current_SNP_name = shared.SNP_ref[min_pC_index];
            LOGGER.i(0, "First SNP", current_SNP_name);
            LOGGER << "p-value: " << scientific << min_pC << fixed << endl;
            LOGGER << "--------------------------------" << endl;

            for (int n : current_list)
                cohorts[n].sumstat_candidate = cohorts[n].sumstat.row(min_pC_index);

            accept_SNP_as_candidate(min_pC_index);
            for (int n : current_list) 
                cohorts[n].save_temp_model();

            iter_num++;
            continue;
        }

        for (int n : current_list)
            cohorts[n].calc_conditional_effects();

        inverse_var_meta_conditional();

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
            
            // temporarily update sumstat_candidate for calculating joint effects
            for (int n : current_list)
                append_row(cohorts[n].sumstat_candidate, cohorts[n].sumstat.row(min_pC_index));

            // calculate joint effects for each cohort
            for (int n : current_list) {
                auto &c = cohorts[n];
                
                if (c.sumstat(min_pC_index, 6) > 1e10) {
                    // LOGGER.w(1, "skipped, adjusted N too large", current_SNP_name);
                    skip_flag = true;
                    break;
                }

                if (!c.calc_R_inverse_forward(min_pC_index)) {
                    // LOGGER.w(1, "skipped, collinearity check failed during R inverse calculation", current_SNP_name);
                    skip_flag = true;
                    break;
                }
                
                if (!c.calc_joint_effects()) {
                    LOGGER.w(1, "skipped, joint se too small", current_SNP_name);
                    skip_flag = true;
                    break;
                }

                // check R2 increment, which does not need to be done in gcta-COJO
                if (!params.if_gcta_COJO && c.R2 < (1+params.R2_threshold) * c.previous_R2) {
                    LOGGER.w(1, "skipped, R2 increment lower than threshold", current_SNP_name);
                    skip_flag = true;
                    break;
                }
            }
            
            // if anything wrong, revert sumstat_candidate and try the next SNP
            if (skip_flag) {
                for (int n : current_list)
                    remove_row(cohorts[n].sumstat_candidate);

                abs_zC(min_pC_index) = -1;
                continue;
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
                accept_SNP_as_candidate(min_pC_index);
                for (int n : current_list)
                    cohorts[n].save_temp_model();
                break;
            } 
            
            // trivial case for removing current SNP
            if (max_pJ_index < fixed_candidate_SNP_num || max_pJ_index == abs_zJ.rows()-1) {
                LOGGER.i(0, "removed, because pJ exceeds p-value threshold", current_SNP_name);
                for (int n : current_list)
                    remove_row(cohorts[n].sumstat_candidate);

                removed_SNP.push_back(min_pC_index);
                auto iter = find(screened_SNP.begin(), screened_SNP.end(), min_pC_index);
                screened_SNP.erase(iter);

                abs_zC(min_pC_index) = -1;
                continue;
            } 

            // backward selection must happen now, save current model as backup
            candidate_SNP_backup = candidate_SNP;
            removed_SNP_backup = removed_SNP;
            for (int n : current_list)
                cohorts[n].save_state();

            LOGGER.i(0, "New candidate SNP added, start backward selection", current_SNP_name);

            accept_SNP_as_candidate(min_pC_index);
            for (int n : current_list)
                cohorts[n].save_temp_model();

            // backward selection
            while (true) {
                LOGGER.i(0, "Candidate SNP removed", shared.SNP_ref[candidate_SNP[max_pJ_index]]);
                remove_SNP_from_candidate(max_pJ_index);

                if (candidate_SNP.size() == 0) {
                    LOGGER.i(0, "Backward selection successful, no candidate SNPs left");
                    break;
                }

                for (int n : current_list) {
                    auto &c = cohorts[n];
                    // do not need to check colinearity because it is guaranteed in forward step
                    c.calc_R_inverse_backward(max_pJ_index);

                    if (!c.calc_joint_effects()) {
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


int MACOJO::set_read_process_output_options(int argc, char** argv)
{   
    CLI::App app;

    app.description(
        "\nMulti-ancestry conditional and joint analysis of GWAS summary statistics\n"
        "https://github.com/light156/multi-ancestry-COJO\n"
        "\nTIPS:\n"
        "  1. Options are basically the same as original GCTA-COJO, while extended to multiple cohorts. \n"
        "  2. Please make sure the paths are paired in --bfile and --cojo-file.\n"
        "     (e.g., --bfile xx1.sumstat xx2.sumstat xx3.sumstat --cojo-file yy1 yy2 yy3)\n"
        "  3. For using genotype data (--bfile), please provide (.bed .bim .fam) files.\n"
        "  4. For using LD matrix (--LD), please provide at least (.bim .ld) files.\n"
        "     If you want to exclude SNPs with too high MAF difference, please provide (.frq) files as well\n"
        "  5. The folders in the output file path (which is after --out) must exist.\n"
        "     (You can easily create them with 'mkdir -p <folder_path>')"
    );

    app.get_formatter()->column_width(70);
    app.get_formatter()->right_column_width(80);

    // options same as GCTA-COJO
    string group1 = "Original GCTA Options/Flags (https://yanglab.westlake.edu.cn/software/gcta/#COJO)";
    app.add_option("--bfile", bfile_list, "PLINK binary file prefix for each cohort (.bed, .bim, .fam)")->required()->group(group1);
    app.add_option("--cojo-file", cojo_file_list, "GWAS summary statistics file for each cohort")->required()->check(CLI::ExistingFile)->group(group1);
    app.add_option("--out", output_name, "Output file path")->required()->group(group1);
    app.add_option("--cojo-wind", params.window_kb, "SNP position window in Kb (-1 for no windows)")->default_val(10000)->group(group1);
    app.add_option("--cojo-p", params.p, "Significance threshold for SNP selection")->default_val(5e-8)->check(CLI::Range(0.0, 1.0))->group(group1);
    app.add_option("--cojo-collinear", params.collinear, "Colinearity threshold (0-0.999)")->default_val(0.9)->check(CLI::Range(0.0, 0.999))->group(group1);
    app.add_option("--maf", params.maf, "Minor Allele Frequency threshold")->default_val(0.01)->check(CLI::Range(1e-5, 0.5))->group(group1);
    app.add_option("--geno", params.missingness, "Missingness threshold")->default_val(1.0)->check(CLI::Range(0.0, 1.0))->group(group1);
    app.add_option("--diff-freq", params.diff_freq, "Frequency difference threshold between sumstat and PLINK (0-1)")->default_val(0.2)->check(CLI::PositiveNumber)->group(group1);
    
    auto *extract_option = app.add_option("--extract", extract_file, "File path of user-given SNPs")->check(CLI::ExistingFile)->group(group1);
    auto *cojo_joint_flag = app.add_flag("--cojo-joint", params.if_cojo_joint, "Only calculate joint effect for provided SNPs and exit")->group(group1);
    cojo_joint_flag->needs(extract_option);

    // multi-cohort specific options
    string group2 = "Manc-COJO Options/Flags";
    auto *fixedSNP_option = app.add_option("--fixed", fixedSNP_file, "File path of fixed candidate SNPs")->check(CLI::ExistingFile)->group(group2);
    auto *R2_option = app.add_option("--R2", params.R2_threshold, "R2 incremental threshold (-1 for no threshold)")->default_val(-1)->group(group2);
    auto *R2back_option = app.add_option("--R2back", params.R2back_threshold, "R2 threshold for backward selection (-1 for no threshold)")->default_val(-1)->group(group2);
    app.add_flag("--freq-mode-and", params.if_freq_mode_and, "Use AND mode for frequency threshold")->group(group2);
    app.add_flag("--MDISA", params.if_MDISA, "Run MDISA after MACOJO")->group(group2);

    // Main logic options
    string group3 = "Main logic Options/Flags";
    int thread_num;
    auto *gcta_flag = app.add_flag("--gcta", params.if_gcta_COJO, "Use GCTA-COJO model selection criteria")->group(group3);
    auto *LD_flag = app.add_flag("--LD", params.if_LD_mode, "Use LD matrix files for analysis")->group(group3);
    auto *remove_NA_flag = app.add_flag("--remove-NA", params.if_remove_NA, "Do not fill NA with mean values in PLINK.bed files")->group(group3);
    auto *iter_option = app.add_option("--iter", params.max_iter_num, "Total iteration number")->default_val(10000)->check(CLI::PositiveNumber)->group(group3);
    app.add_option("--thread-num", thread_num, "Number of threads to use")->default_val(1)->check(CLI::Range(1,10))->group(group3);

    remove_NA_flag->excludes(gcta_flag);
    remove_NA_flag->excludes(LD_flag);

    for (auto *opt : app.get_options())
        opt->multi_option_policy(CLI::MultiOptionPolicy::Throw);

    app.callback([&]() {
        if (bfile_list.size() != cojo_file_list.size())
            throw CLI::ValidationError("Filepath numbers after --bfile and --cojo-file must be the same, please check");

        if (bfile_list.size() == 1 && (params.if_freq_mode_and || params.if_MDISA))
            throw CLI::ValidationError("--freq-mode-and/--MDISA are invalid options with single cohort");

        if (*cojo_joint_flag && (*fixedSNP_option || *iter_option || *R2_option || *R2back_option))
            throw CLI::ValidationError("--cojo-joint cannot be used with --fixed/--iter/--R2/--R2back, because it makes no sense");
        
        if (*cojo_joint_flag && bfile_list.size() > 1)
            throw CLI::ValidationError("--cojo-joint should only be used with single cohort");

        if (*gcta_flag && (*R2_option || *R2back_option))
            throw CLI::ValidationError("--gcta cannot be used with --R2/--R2back, because GCTA-COJO does not check R2 increment");
    });

    try { 
        app.parse(argc, argv); 
    } catch (const CLI::CallForHelp &e) {
        LOGGER << app.help() << endl;
        return -1;
    } catch (const CLI::ParseError &e) {
        LOGGER << app.help() << endl;
        LOGGER.e(0, e.what());
    }

    omp_set_num_threads(thread_num);
    params.window_size = (params.window_kb < 0 ? INT_MAX : params.window_kb * 1000);
    params.iter_collinear_threshold = 1.0 / (1.0 - params.collinear);
    LOGGER.open(output_name + ".log");

    for (size_t i = 0; i < bfile_list.size(); i++)
        cohorts.emplace_back(params, shared);

    if (!fixedSNP_file.empty()) {
        read_SNP_only(fixedSNP_file, fixed_candidate_SNP);
        LOGGER.i(0, "fixed candidate SNPs provided by the user for analysis", to_string(fixed_candidate_SNP.size()));
    }

    // output_user_hyperparameters
    LOGGER << "=========== MACOJO CONFIGURATION ===========" << endl
            << (params.if_cojo_joint ? "COJO-Joint analysis only\n" : "COJO iterative selection\n")
            << "p-value Threshold: " << params.p << endl
            << "Collinearity threshold: " << params.collinear << endl
            << "SNP position window (+/-): " << params.window_kb << "Kb" << endl
            << "Sumstat Minor Allele Frequency threshold: " << params.maf << endl
            << "Genotype Missingness threshold: " << params.missingness << endl
            << "Frequency difference threshold between sumstat and PLINK: " << params.diff_freq << endl
            << (*R2_option ? "R2 incremental threshold: " + to_string(params.R2_threshold) + "\n" : "")
            << (*R2back_option ? "R2 incremental threshold (backward): " + to_string(params.R2back_threshold) + "\n" : "")
            << (cohorts.size() > 1 ? 
                (params.if_freq_mode_and ? "SNP frequency mode: AND\n" : "SNP frequency mode: OR\n") : "")
            << (params.if_MDISA ? "Run MDISA after MACOJO\n" : "")
            << (params.if_LD_mode ? "Use .ld files for calculation\n" : "")
            << (params.if_remove_NA ? "Do not fill NA with mean genotype values\n" : "")
            << (params.if_gcta_COJO ? "Use GCTA-COJO model selection criteria\n" : "")
            << "===========================================" << endl << endl;

    return 0;
}
