#include "include/macojo.h"


// Map user-provided SNP names to internal indices.
void MACOJO::initialize_candidate_SNP(const vector<string>& given_list)
{
    for (const auto& SNP_name : given_list) {
        auto iter = fast_lookup(shared.goodSNP_table, SNP_name);
        if (iter == shared.goodSNP_table.end()){
            // LOGGER.w("not found in current SNP list and will be ignored", SNP_name);
            continue;
        }
        candidate_SNP.push_back(iter->second);
    }
}


// Build LD matrix (R) and its inverse from scratch for a candidate list.
int Cohort::calc_R_inv_from_SNP_list(const vector<int>& SNP_list, string mode) 
{   
    if (SNP_list.size() == 0) 
        LOGGER.e("No candidate SNPs for calculation, program exit!");

    int total_num = SNP_list.size();
    R_post = MatrixXd::Identity(total_num, total_num);

    for (int i = 0; i < total_num; i++) {
        for (int j = 0; j < i; j++) {
            if (abs(shared.SNP_pos_ref[SNP_list[i]] - shared.SNP_pos_ref[SNP_list[j]]) < params.window_size) {
                if (params.if_LD_mode)
                    R_post(i, j) = LD_matrix(SNP_list[i], SNP_list[j]);
                else
                    R_post(i, j) = genotype.calc_inner_product(SNP_list[i], SNP_list[j], mode == "removeNA");

                R_post(j, i) = R_post(i, j);
            }
        }
    }
    
    R_inv_post = R_post.ldlt().solve(MatrixXd::Identity(total_num, total_num));

    if (mode == "GCTA") {
        MatrixXd R_post_gcta = MatrixXd::Zero(total_num, total_num);

        for (int i = 0; i < total_num; i++) {
            R_post_gcta(i, i) = sumstat(SNP_list[i], 6);

            for (int j = 0; j < i; j++) {
                R_post_gcta(i, j) = R_post(i, j) * min(sumstat(SNP_list[i], 4), sumstat(SNP_list[j], 4)) * \
                    sqrt(sumstat(SNP_list[i], 5) * sumstat(SNP_list[j], 5));
                R_post_gcta(j, i) = R_post_gcta(i, j);
            }
        }

        R_inv_post_gcta = R_post_gcta.ldlt().solve(MatrixXd::Identity(total_num, total_num));
    }
    
    int delete_index;

    if (R_inv_post.cwiseAbs().diagonal().maxCoeff(&delete_index) > params.iter_collinear_threshold)
        return delete_index;

    if (calc_joint_effects(SNP_list, mode) == -1) {
        beta_var.minCoeff(&delete_index);
        return delete_index;
    }
    
    save_temp_model();
    return -1;
}


// It is called "check", but it builds R and R^{-1} from scratch and removes
// collinear SNPs until the set is stable.
bool MACOJO::check_candidate_SNP_collinearity(string mode)
{
    bool success_flag = true;

    // get R_inv from scratch
    while (true) {
        int i;
        for (i = 0; i < current_list.size(); i++) {
            int delete_index = cohorts[current_list[i]].calc_R_inv_from_SNP_list(candidate_SNP, mode);
            if (delete_index == -1) continue;

            LOGGER.w("removed due to collinearity", shared.goodSNP_table[candidate_SNP[delete_index]].first);
            collinear_SNP.push_back(candidate_SNP[delete_index]);
            candidate_SNP.erase(candidate_SNP.begin() + delete_index);
            success_flag = false;
            break;
        }
        if (i == current_list.size()) break;
    }

    return success_flag;
}


// Main dispatcher for joint/conditional/stepwise selection workflows.
void MACOJO::entry_function()
{   
    LOGGER.i("Calculation started!");

    for (int n = 0; n < cohorts.size(); n++)
        current_list.push_back(n);

    if (params.if_joint_mode) {
        for (int i = 0; i < shared.goodSNP_table.size(); i++) {
            if (shared.goodSNP_table[i].second != -1)
                candidate_SNP.push_back(i);
        }

        output_jma(params.output_name);
        LOGGER.i("Calculation finished for joint analysis!");
        return;
    }  

    if (params.if_cond_mode) {
        initialize_candidate_SNP(params.cojo_cond_SNPs);
        LOGGER.i("valid conditional SNPs provided by the user", candidate_SNP.size()); 
        if (candidate_SNP.size() == 0) return;

        if (!check_candidate_SNP_collinearity(params.effect_size_mode))
            LOGGER.i("valid conditional SNPs after checking collinearity", candidate_SNP.size());

        output_cma(params.output_name);
        LOGGER.i("Calculation finished for conditional analysis!");
        return;
    }
    
    // COJO stepwise iterative selection, first deal with fixed candidate SNPs if provided
    if (!params.fix_SNPs.empty()) {
        initialize_candidate_SNP(params.fix_SNPs);
        if (candidate_SNP.size() > 0) {
            LOGGER.i("valid fixed SNPs provided by the user", candidate_SNP.size());

            if (!check_candidate_SNP_collinearity(params.slct_mode))
                LOGGER.i("valid fixed SNPs after checking collinearity", candidate_SNP.size());

            for (int idx : candidate_SNP)
                active_mask[idx] = 0;

            // initialize r and r_gcta from scratch
            for (auto &c : cohorts) {
                for (int index : candidate_SNP)
                    c.append_r(active_mask, index, params.slct_mode);
            }
        } else
            LOGGER.i("No valid fixed SNPs found, proceed with stepwise selection from zero");
    }

    fixed_candidate_SNP_num = candidate_SNP.size();
    slct_loop();
    output_jma(params.output_name);

    if (cohorts.size() > 1 && params.if_MDISA) {
        current_list.resize(1);

        // backup SNPs for MDISA, Manc-COJO candidate SNPs will be fixed
        fixed_candidate_SNP_num = candidate_SNP.size();
        candidate_SNP_backup = candidate_SNP;

        for (size_t n = 0; n < cohorts.size(); n++) {
            LOGGER << endl << "MDISA for Cohort " << n + 1 << endl;
            current_list[0] = n;
                
            candidate_SNP = candidate_SNP_backup;
            vector<int>().swap(collinear_SNP);
            vector<int>().swap(backward_SNP);

            active_mask.assign(shared.goodSNP_table.size(), 1);
            for (int idx : bad_SNP)
                active_mask[idx] = 0;
            for (int idx : candidate_SNP)
                active_mask[idx] = 0;
                
            // calculate r and r_gcta from scratch
            cohorts[n].r.resize(0, 0);
            cohorts[n].r_gcta.resize(0, 0);
            for (int index : candidate_SNP)
                cohorts[n].append_r(active_mask, index, params.slct_mode);

            // calculte R from scratch
            cohorts[n].calc_R_inv_from_SNP_list(candidate_SNP, params.slct_mode);

            slct_loop();
            output_jma(params.output_name + ".MDISA.cohort" + to_string(n + 1));
        }
    }

    LOGGER.i("All calculations finished!");
}


// Write conditional results for all active SNPs (excluding masked SNPs).
// Must happen after checking candidate SNP collinearity.
void MACOJO::output_cma(string savename) 
{   
    active_mask.assign(shared.goodSNP_table.size(), 1);
    for (int idx : bad_SNP)
        active_mask[idx] = 0;
    for (int idx : candidate_SNP)
        active_mask[idx] = 0;
    for (int idx : collinear_SNP)
        active_mask[idx] = 0;
    for (int idx : backward_SNP)
        active_mask[idx] = 0;

    size_t active_count = 0;
    for (char v : active_mask)
        active_count += (v != 0);
    if (active_count == 0) {
        LOGGER.i("No screened SNPs, conditional analysis output file will not be generated");
        return;
    }

    // calculate r and r_gcta from scratch              
    if (params.if_cond_mode || params.slct_mode != params.effect_size_mode) {
        for (int n : current_list) {
            cohorts[n].r.resize(0, 0);
            cohorts[n].r_gcta.resize(0, 0);
            for (int index : candidate_SNP)
                cohorts[n].append_r(active_mask, index, params.effect_size_mode);
        }
    }

    // since r and R_inv_pre are ready, calculate conditional effects
    for (int n : current_list) 
        cohorts[n].calc_cond_effects(candidate_SNP, params.effect_size_mode);
    
    vector<pair<int, int>> SNP_ref_order_pair;
    SNP_ref_order_pair.reserve(shared.bp_order.size());
    for (int idx : shared.bp_order) {
        if (active_mask[idx])
            SNP_ref_order_pair.emplace_back(shared.SNP_pos_ref[idx], idx);
    }

    output_inverse_var_meta(savename, SNP_ref_order_pair, false);
}


// Write joint results for current candidate SNPs.
void MACOJO::output_jma(string savename) 
{   
    if (candidate_SNP.size() == 0) {
        LOGGER.i("No candidate SNPs, joint analysis output file will not be generated");
        return;
    }

    // calculate from scratch              
    if (params.if_joint_mode || params.slct_mode != params.effect_size_mode) {
        if (!check_candidate_SNP_collinearity(params.effect_size_mode))
            LOGGER.i("valid joint SNPs after checking collinearity", candidate_SNP.size());
    } 

    // we may need output R, so a duplicate calculation simplifies the code structure here
    for (int n : current_list) 
        cohorts[n].calc_R_inv_from_SNP_list(candidate_SNP, params.effect_size_mode);

    vector<int> cand_pos(shared.goodSNP_table.size(), -1);
    for (int i = 0; i < candidate_SNP.size(); i++)
        cand_pos[candidate_SNP[i]] = i;

    vector<pair<int, int>> SNP_ref_order_pair;
    SNP_ref_order_pair.reserve(candidate_SNP.size());
    for (int idx : shared.bp_order) {
        int pos = cand_pos[idx];
        if (pos != -1)
            SNP_ref_order_pair.emplace_back(shared.SNP_pos_ref[idx], pos);
    }

    output_inverse_var_meta(savename, SNP_ref_order_pair, true);

    if (params.if_output_all) 
        output_ld_matrix(savename, SNP_ref_order_pair);

    if (params.if_output_all && !params.if_joint_mode) 
        output_cma(savename);
}


void MACOJO::output_inverse_var_meta(string savename, const vector<pair<int, int>>& SNP_ref_order_pair, bool if_joint) 
{   
    savename += if_joint ? ".jma.cojo" : ".cma.cojo";
    char mode = if_joint ? 'J' : 'C';

    bool if_exist_file = false;
    bool is_first_chr = (params.curr_chr == params.chr_list.front());

    ifstream check(savename.c_str());
    if (check) if_exist_file = true;
    check.close();

    ofstream maCOJO;

    if (if_exist_file && !is_first_chr)
        maCOJO.open(savename.c_str(), ios::app);
    else 
        maCOJO.open(savename.c_str());
    
    if (!maCOJO) LOGGER.e("Cannot open the file [" + savename + "] to write");
    maCOJO.precision(10);

    if (!if_exist_file || is_first_chr) {
        // write header
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
    }

    ArrayXd bma, se2ma, abs_zma;
    if (current_list.size() > 1)
        inverse_var_meta(bma, se2ma, abs_zma);

    for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) {
        
        int info_index = if_joint ? candidate_SNP[iter->second] : iter->second;
        maCOJO << params.curr_chr
                << "\t" << shared.goodSNP_table[info_index].first
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


// Write the LD matrix for candidate SNPs (if requested).
void MACOJO::output_ld_matrix(string savename, const vector<pair<int, int>>& SNP_ref_order_pair) 
{   
    savename += ".ldr.cojo";

    bool if_exist_file = false;
    bool is_first_chr = (params.curr_chr == params.chr_list.front());

    ifstream check(savename.c_str());
    if (check) if_exist_file = true;
    check.close();

    ofstream ldrCOJO;

    if (if_exist_file && !is_first_chr)
        ldrCOJO.open(savename.c_str(), ios::app);
    else 
        ldrCOJO.open(savename.c_str());

    if (!ldrCOJO) LOGGER.e("Cannot open the file [" + savename + "] to write");
    ldrCOJO.precision(10);
    
    // calc_R_inv_from_SNP_list has just been called so R_post is ready
    for (int n : current_list) {
        if (params.chr_list.size() > 1)
            ldrCOJO << "# Chromosome " << params.curr_chr << "\n";
        
        if (current_list.size() > 1)
            ldrCOJO << "# Cohort " << n + 1 << "\n";
        
        ldrCOJO << "SNP";

        for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) 
            ldrCOJO << "\t" << shared.goodSNP_table[candidate_SNP[iter->second]].first;

        ldrCOJO << "\n";

        auto &c = cohorts[n];
        int total_num = c.R_post.rows();

        for (auto iter1 = SNP_ref_order_pair.begin(); iter1 != SNP_ref_order_pair.end(); iter1++) {
            ldrCOJO << shared.goodSNP_table[candidate_SNP[iter1->second]].first;
            for (auto iter2 = SNP_ref_order_pair.begin(); iter2 != SNP_ref_order_pair.end(); iter2++)
                ldrCOJO << "\t" << c.R_post(iter1->second, iter2->second);

            ldrCOJO << "\n";
        }
    }

    ldrCOJO.close();
    LOGGER.i("Results saved into [" + savename + "]");
}


// Write SNPs filtered out for quality control reasons.
void MACOJO::output_bad_SNP(string savename) 
{   
    if (shared.bad_SNP_dict.size() == 0) {
        LOGGER.i("No bad SNPs detected");
        return;
    }

    savename += ".badsnps";

    bool if_exist_file = false;
    bool is_first_chr = (params.curr_chr == params.chr_list.front());

    ifstream check(savename.c_str());
    if (check) if_exist_file = true;
    check.close();

    ofstream badSNPout;

    if (if_exist_file && !is_first_chr)
        badSNPout.open(savename.c_str(), ios::app);
    else 
        badSNPout.open(savename.c_str());

    if (!badSNPout) LOGGER.e("Cannot open the file [" + savename + "] to write");

    badSNPout << "# Chromosome " << params.curr_chr << "\n";

    sort(shared.bad_SNP_dict.begin(), shared.bad_SNP_dict.end());
    shared.bad_SNP_dict.erase(unique(shared.bad_SNP_dict.begin(), shared.bad_SNP_dict.end()), shared.bad_SNP_dict.end());

    BadSnpReason current = shared.bad_SNP_dict.begin()->first;
    badSNPout << "[" << reason_to_string(current) << "]\n";

    for (auto iter = shared.bad_SNP_dict.begin(); iter != shared.bad_SNP_dict.end(); iter++) {
        if (iter->first != current) {
            current = iter->first;
            badSNPout << "\n[" << reason_to_string(current) << "]\n";
        }
        badSNPout << "  " << iter->second;
    }

    badSNPout << "\n";
    badSNPout.close();
    LOGGER.i("List of bad SNPs saved into [" + savename + "]");
}
