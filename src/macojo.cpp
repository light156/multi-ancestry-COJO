#include "macojo.h"


void MACOJO::initialize_candidate_SNP(string filename)
{
    vector<string> temp;
    skim_file(filename, temp, false, true, false);

    for (const auto& SNP_name : temp) {
        auto iter = shared.goodSNP_index_map.find(SNP_name);
        if (iter == shared.goodSNP_index_map.end()){
            LOGGER.w("not found in + [ " + filename + " ] and will be ignored", SNP_name);
            continue;
        }
        candidate_SNP.push_back(iter->second);
        screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), iter->second));
    }

    if (candidate_SNP.size() == 0)
        LOGGER.e("None of the provided SNPs exist in data");
}


void MACOJO::entry_function()
{   
    LOGGER.i("Calculation started!");

    for (size_t n = 0; n < cohorts.size(); n++)
        current_list.push_back(n);

    if (params.if_joint_mode) {
        candidate_SNP = screened_SNP; // screened SNP cannot be empty here, checked during reading files
        output_jma(params.output_name);
        LOGGER.i("Calculation finished for joint analysis, program exit!");
        return;
    }  

    if (params.if_cond_mode) {
        initialize_candidate_SNP(params.cond_file);
        LOGGER.i("effective conditional SNPs provided by the user", to_string(candidate_SNP.size()));

        output_cma(params.output_name);
        LOGGER.i("Calculation finished for conditional analysis, program exit!");
        return;
    }
    
    // COJO stepwise iterative selection, first deal with fixed candidate SNPs if provided
    if (!params.fixedSNP_file.empty()) {
        initialize_candidate_SNP(params.fixedSNP_file);
        LOGGER.i("effective fixed candidate SNPs provided by the user", to_string(fixed_candidate_SNP_num));

        bool success_flag = true;

        for (auto &c : cohorts) {
            for (int index : candidate_SNP)
                c.append_r(screened_SNP, index, params.slct_mode);

            success_flag &= c.calc_R_inv_from_SNP_list(candidate_SNP, params.slct_mode);
            c.save_temp_model();

            if (candidate_SNP.size() == 0) {
                LOGGER.i("No valid fixed candidate SNPs after checking collinearity, program exit!");
                return;
            }
        }

        if (!success_flag) {
            LOGGER.i("Proceeding to stepwise selection with remaining fixed candidate SNPs");
            params.output_name += ".warning";
        }
    }

    fixed_candidate_SNP_num = candidate_SNP.size();
    slct_loop();
    output_jma(params.output_name);
    if (params.if_output_all)
        output_cma(params.output_name);

    if (cohorts.size() > 1 && params.if_MDISA) {
        current_list.resize(1);

        // backup SNPs for MDISA, Manc-COJO candidate SNPs will be fixed
        fixed_candidate_SNP_num = candidate_SNP.size();
        candidate_SNP_backup = candidate_SNP;

        for (size_t n = 0; n < cohorts.size(); n++) {
            LOGGER << "MDISA for Cohort " << n + 1 << endl;
            current_list[0] = n;
            
            candidate_SNP = candidate_SNP_backup;
            screened_SNP.insert(screened_SNP.end(), collinear_SNP.begin(), collinear_SNP.end());
            screened_SNP.insert(screened_SNP.end(), backward_SNP.begin(), backward_SNP.end());
            vector<int>().swap(collinear_SNP);
            vector<int>().swap(backward_SNP);

            slct_loop();
            output_jma(params.output_name + ".MDISA.cohort" + to_string(n + 1));
            if (params.if_output_all)
                output_cma(params.output_name + ".MDISA.cohort" + to_string(n + 1));
        }
    }

    LOGGER.i("All calculations finished!");
}


void MACOJO::output_cma(string savename) 
{   
    if (screened_SNP.size() == 0) {
        LOGGER.i("No screened SNPs, conditional analysis output file will not be generated");
        return;
    }

    // calculate from scratch
    if (params.if_cond_mode || params.slct_mode != params.effect_size_mode) {
        bool success_flag = true;

        for (int n : current_list) {
            auto& c = cohorts[n];

            // get r and r_gcta
            for (int index : candidate_SNP)
                c.append_r(screened_SNP, index, params.effect_size_mode);
            
            // get R_inv_pre and R_inv_pre_gcta from R_inv_post and R_inv_post_gcta
            success_flag &= c.calc_R_inv_from_SNP_list(candidate_SNP, params.effect_size_mode);
            c.save_temp_model();

            if (candidate_SNP.size() == 0) {
                LOGGER.i("No valid candidate SNPs after checking collinearity, conditional analysis output file will not be generated");
                return;
            }
        }

        if (!success_flag) {
            LOGGER.i("Some candidate SNPs were removed, please be aware of this");
            savename += ".warning";
        }
    }

    // since r and R_inv_pre are ready, calculate conditional effects
    for (int n : current_list) 
        cohorts[n].calc_cond_effects(candidate_SNP, params.effect_size_mode);

    inverse_var_meta(bC, se2C, abs_zC);
    
    map<int, int> SNP_ref_order_pair;
    for (int index : screened_SNP)
        SNP_ref_order_pair.insert(make_pair(shared.SNP_pos_ref[index], index));

    output_inverse_var_meta(savename + ".cma.cojo", 'C', SNP_ref_order_pair, bC, se2C, abs_zC);
}


void MACOJO::output_jma(string savename) 
{   
    if (candidate_SNP.size() == 0) {
        LOGGER.i("No candidate SNPs, joint analysis output file will not be generated");
        return;
    }

    // calculate from scratch for simplicity
    bool success_flag = true;

    for (int n : current_list) {
        auto& c = cohorts[n];
        
        success_flag &= c.calc_R_inv_from_SNP_list(candidate_SNP, params.effect_size_mode);
        c.save_temp_model();

        if (candidate_SNP.size() == 0) {
            LOGGER.i("No valid candidate SNPs after checking collinearity, joint analysis output file will not be generated");
            return;
        }
    }

    if (!success_flag) {
        LOGGER.i("Some candidate SNPs were removed, please be aware of this");
        savename += ".warning";
    }

    inverse_var_meta(bJ, se2J, abs_zJ);

    map<int, int> SNP_ref_order_pair;
    for (int num = 0; num < candidate_SNP.size(); num++)
        SNP_ref_order_pair.insert(make_pair(shared.SNP_pos_ref[candidate_SNP[num]], num));

    output_inverse_var_meta(savename + ".jma.cojo", 'J', SNP_ref_order_pair, bJ, se2J, abs_zJ);

    if (params.if_output_all)  {
        vector<int> ordered_candidate;
        for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++)
            ordered_candidate.push_back(iter->second); 

        if (current_list.size() == 1)
            output_ld_matrix(savename + ".ldr.cojo", ordered_candidate, cohorts[current_list[0]]);
        else
            for (int n : current_list)
                output_ld_matrix(savename + ".cohort" + to_string(n+1) + ".ldr.cojo", ordered_candidate, cohorts[n]);
    }   
}


void MACOJO::output_inverse_var_meta(string savename, char mode, const map<int, int>& SNP_ref_order_pair, 
                                const ArrayXd& bma, const ArrayXd& se2ma, const ArrayXd& abs_zma) 
{
    ofstream maCOJO(savename.c_str());
    if (!maCOJO) LOGGER.e("Cannot open the file [" + savename + "] to write");

    maCOJO.precision(12);
    maCOJO << "Chr\tSNP\tbp\tA1\tA2";

    if (current_list.size() == 1)
        maCOJO << "\tfreq"
                << "\tb"
                << "\tse"
                << "\tp"
                << "\tn"
                << "\tb" << mode 
                << "\tb" << mode <<"_se"
                << "\tp" << mode;
    else {
        for (int n : current_list) {
            maCOJO << "\tfreq." << n+1 
                    << "\tb." << n+1 
                    << "\tse." << n+1 
                    << "\tp." << n+1 
                    << "\tn." << n+1 
                    << "\tb" << mode << "." << n+1 
                    << "\tb" << mode << "_se." << n+1 
                    << "\tp" << mode << "." << n+1;
        }
        
        maCOJO << "\tb" << mode <<".ma"  
                << "\tb" << mode << "_se.ma" 
                << "\tp" << mode << ".ma";
    }

    maCOJO << "\n";

    for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) {
        
        int info_index = (mode == 'C') ? iter->second : candidate_SNP[iter->second];
        maCOJO << shared.chr_ref[info_index]
                << "\t" << shared.SNP_ref[info_index]
                << "\t" << iter->first
                << "\t" << shared.A1_ref[info_index]
                << "\t" << shared.A2_ref[info_index];

        for (int n : current_list) {
            auto &c = cohorts[n];
            
            // cols 0:b, 1:se2, 2:p, 3:freq, 4:n, 5: b, 6:se, 7:p
            maCOJO << "\t" << c.sumstat(info_index, 3) 
                    << "\t" << c.sumstat(info_index, 0)
                    << "\t" << sqrt(c.sumstat(info_index, 1))
                    << "\t" << c.sumstat(info_index, 2)
                    << "\t" << c.sumstat(info_index, 4)
                    << "\t" << c.beta(iter->second)
                    << "\t" << sqrt(c.beta_var(iter->second))
                    << "\t" << erfc(abs(c.beta(iter->second)) / sqrt(c.beta_var(iter->second)) / sqrt(2));
        }

        if (current_list.size() > 1) {
            maCOJO << "\t" << bma(iter->second) 
                    << "\t" << sqrt(se2ma(iter->second))
                    << "\t" << erfc(abs_zma(iter->second) / sqrt(2));
        } 
        
        maCOJO << "\n";
    }

    maCOJO.close();
    LOGGER.i("Results saved into [" + savename + "]");
}


void MACOJO::output_ld_matrix(string savename, const vector<int>& ordered_candidate, const Cohort& c) 
{   
    ofstream ldrCOJO(savename);
    if (!ldrCOJO) LOGGER.e("Cannot open the file [" + savename + "] to write");

    ldrCOJO.precision(12);
    ldrCOJO << "SNP";

    for (int index : ordered_candidate) 
        ldrCOJO << "\t" << shared.SNP_ref[candidate_SNP[index]];

    ldrCOJO << "\n";

    // calc_R_inv_from_SNP_list has just been called so R_post is ready
    int total_num = c.R_post.rows();

    for (int index_i : ordered_candidate) {
        ldrCOJO << shared.SNP_ref[candidate_SNP[index_i]];
        for (int index_j : ordered_candidate)
            ldrCOJO << "\t" << c.R_post(index_i, index_j);

        ldrCOJO << "\n";
    }

    ldrCOJO.close();
    LOGGER.i("Results saved into [" + savename + "]");
}
