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

    for (int idx : candidate_SNP)
        active_mask[idx] = 0;
}


// Build LD matrix (R) and its inverse from scratch for a candidate list.
int Cohort::calc_R_inv_and_joint(const vector<int>& SNP_list, int fixed_candidate_num, string mode)
{
    // Empty list: reset R_inv_pre so that calc_R_inv_forward treats the next
    // SNP as the very first candidate (base case: R_inv_pre.rows() == 0).
    // This can happen when --fix-drop removes all fixed SNPs.
    if (SNP_list.empty()) {
        R_inv_pre.resize(0, 0);
        R_inv_pre_gcta.resize(0, 0);
        previous_R2 = 0;
        return -1;
    }

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
    double max_abs_diag;

    const int fixed_dim = min(max(fixed_candidate_num, 0), static_cast<int>(R_inv_post.rows()));
    const int free_dim = static_cast<int>(R_inv_post.rows()) - fixed_dim;

    if (fixed_dim > 0) {
        max_abs_diag = R_inv_post.topLeftCorner(fixed_dim, fixed_dim).cwiseAbs().diagonal().maxCoeff(&delete_index);
        if (max_abs_diag > 1.0 / (1.0 - params.fix_collinear))
            return delete_index;
    }

    if (free_dim > 0) {
        max_abs_diag = R_inv_post.bottomRightCorner(free_dim, free_dim).cwiseAbs().diagonal().maxCoeff(&delete_index);
        if (max_abs_diag > 1.0 / (1.0 - params.collinear))
            return delete_index + fixed_dim;
    }

    if (calc_joint_effects(SNP_list, mode) == -1) {
        beta_var.minCoeff(&delete_index);
        return delete_index;
    }

    R_inv_pre = R_inv_post;
    R_inv_pre_gcta = R_inv_post_gcta;
    previous_R2 = R2;
    return -1;
}


// It is called "check", but it builds R and R^{-1} and calculate effect sizes
// from scratch and removes collinear SNPs until the set is stable.
void MACOJO::check_candidate_SNP_collinearity(string mode)
{
    // get R_inv from scratch
    while (true) {
        size_t i;
        for (i = 0; i < cohorts.size(); i++) {
            int delete_index = cohorts[i].calc_R_inv_and_joint(candidate_SNP, fixed_candidate_SNP_num, mode);
            if (delete_index == -1) continue;

            LOGGER.w("removed due to collinearity in Cohort " + to_string(i+1), shared.goodSNP_table[candidate_SNP[delete_index]].first);
            collinear_SNP.push_back(candidate_SNP[delete_index]);
            candidate_SNP.erase(candidate_SNP.begin() + delete_index);
            break;
        }
        if (i == cohorts.size()) break;
    }
}


void MACOJO::prepare_fixed_SNP()
{
    while(true) {
        check_candidate_SNP_collinearity(params.slct_mode);
        if (candidate_SNP.empty()) break;

        int max_pJ_index;
        double max_pJ;

        // calculate joint p-value
        inverse_var_meta(bJ, se2J, abs_zJ);
        max_pJ = erfc(abs_zJ.minCoeff(&max_pJ_index) / sqrt(2));
        LOGGER << "In fixed SNPs, max joint pJ: " << max_pJ << endl;
        if (max_pJ <= params.fix_p_value) break;

        if (!params.if_fix_drop) {
            LOGGER.w("p-value for fixed SNPs exceeds threshold, which may affect the following selection process. "
                "Consider using --fix-drop to iteratively remove SNPs, or make sure this is what you want");
            break;
        }

        LOGGER.i("removed, because user specifies --fix-drop" , shared.goodSNP_table[candidate_SNP[max_pJ_index]].first);
        bad_SNP.push_back(candidate_SNP[max_pJ_index]);
        candidate_SNP.erase(candidate_SNP.begin() + max_pJ_index);
    }

    // Build r for surviving fixed SNPs so slct_loop() starts with the correct
    // LD columns. If candidate_SNP is empty (all fixed SNPs dropped by --fix-drop),
    // this loop is a no-op and R_inv_pre stays 0x0 (already reset by
    // check_candidate_SNP_collinearity -> calc_R_inv_and_joint when given empty list).
    for (auto &c : cohorts)
        for (int index : candidate_SNP)
            c.append_r(active_mask, index, params.slct_mode);
}

// Main dispatcher for joint/conditional/stepwise selection workflows.
void MACOJO::entry_function()
{
    LOGGER.i("Calculation started!");

    if (params.if_joint_mode) {
        for (size_t i = 0; i < shared.goodSNP_table.size(); i++) {
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

        output_cma(params.output_name, true);
        LOGGER.i("Calculation finished for conditional analysis!");
        return;
    }

    // COJO stepwise iterative selection, first deal with fixed candidate SNPs if provided
    if (!params.fix_SNPs.empty()) {
        initialize_candidate_SNP(params.fix_SNPs);
        LOGGER.i("valid fixed SNPs provided by the user", candidate_SNP.size());

        LOGGER << endl << "Fixed p-value threshold: " << params.fix_p_value << endl;
        LOGGER << "Collinearity threshold for fixed SNPs: " << params.fix_collinear << endl;
        LOGGER << "Collinearity threshold between fixed SNPs and selected SNPs: " << params.fix_cojo_collinear << endl << endl;

        fixed_candidate_SNP_num = candidate_SNP.size();
        prepare_fixed_SNP();
        LOGGER.i("fixed SNPs retained for stepwise selection", candidate_SNP.size());
    }

    fixed_candidate_SNP_num = candidate_SNP.size();
    slct_loop();
    output_jma(params.output_name);
    LOGGER.i("All calculations finished!");
}


// Write conditional results for all active SNPs (excluding masked SNPs).
// Must happen after checking candidate SNP collinearity.
void MACOJO::output_cma(string savename, bool force_rebuild_r)
{
    if (params.if_cond_mode) {
        check_candidate_SNP_collinearity(params.effect_size_mode);
        LOGGER.i("valid conditional SNPs after checking collinearity", candidate_SNP.size());
    }

    if (candidate_SNP.empty()) {
        LOGGER.i("No conditional SNPs, conditional analysis output file will not be generated");
        return;
    }

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
    if (force_rebuild_r) {
        for (auto& c : cohorts) {
            c.r.resize(0, 0);
            c.r_gcta.resize(0, 0);
            for (int index : candidate_SNP)
                c.append_r(active_mask, index, params.effect_size_mode);
        }
    }

    // since r and R_inv_pre are ready, calculate conditional effects
    for (auto& c : cohorts)
        c.calc_cond_effects(candidate_SNP, params.effect_size_mode);

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
    size_t cand_before = candidate_SNP.size();

    // check_candidate_SNP_collinearity will calculate R and joint effects for all candidate SNPs from scratch
    // we may need output R for ld matrix, or handle with case when the last iteration try to select some SNPs
    // but did not pass check (this causes a bug in GCTA), so calculating from scratch is the best strategy here
    check_candidate_SNP_collinearity(params.effect_size_mode);
    LOGGER.i("valid joint SNPs after checking collinearity", candidate_SNP.size());

    if (candidate_SNP.empty()) {
        LOGGER.i("No joint SNPs, joint analysis output file will not be generated");
        return;
    }

    vector<int> cand_pos(shared.goodSNP_table.size(), -1);
    for (size_t i = 0; i < candidate_SNP.size(); i++)
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

    if (params.if_output_all && !params.if_joint_mode) {
        bool force_rebuild_r = (
            (candidate_SNP.size() != cand_before) ||
            (params.slct_mode != params.effect_size_mode)
        );
        output_cma(savename, force_rebuild_r);
    }
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

        if (cohorts.size() == 1)
            maCOJO << "\tfreq"
                    << "\tb"
                    << "\tse"
                    << "\tp"
                    << "\tn"
                    << "\tb" << mode
                    << "\tb" << mode <<"_se"
                    << "\tp" << mode;
        else {
            for (size_t n = 0; n < cohorts.size(); n++) {
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

            if (params.if_hetero_report)
                maCOJO << "\tCochran_Q";
        }

        maCOJO << "\n";
    }

    ArrayXd bma, se2ma, abs_zma;
    if (cohorts.size() > 1)
        inverse_var_meta(bma, se2ma, abs_zma);

    for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) {

        int info_index = if_joint ? candidate_SNP[iter->second] : iter->second;
        maCOJO << params.curr_chr
                << "\t" << shared.goodSNP_table[info_index].first
                << "\t" << iter->first
                << "\t" << shared.A1_ref[info_index]
                << "\t" << shared.A2_ref[info_index];

        for (auto& c : cohorts) {
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

        if (cohorts.size() > 1) {
            maCOJO << "\t" << bma(iter->second)
                    << "\t" << sqrt(se2ma(iter->second))
                    << "\t" << erfc(abs_zma(iter->second) / sqrt(2));

            if (params.if_hetero_report) {
                // Cochran's Q = sum_i [ (b_i - b_meta)^2 / se2_i ]
                double Q = 0.0;
                for (auto& c : cohorts) {
                    double diff = c.beta(iter->second) - bma(iter->second);
                    Q += diff * diff / c.beta_var(iter->second);
                }
                maCOJO << "\t" << Q;
            }
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

    // calc_R_inv_and_joint has just been called so R_post is ready
    for (const auto& c : cohorts) {
        if (params.chr_list.size() > 1)
            ldrCOJO << "# Chromosome " << params.curr_chr << "\n";

        if (cohorts.size() > 1)
            ldrCOJO << "# Cohort " << c.cohort_index + 1 << "\n";

        ldrCOJO << "SNP";

        for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++)
            ldrCOJO << "\t" << shared.goodSNP_table[candidate_SNP[iter->second]].first;

        ldrCOJO << "\n";

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
    if (shared.bad_SNP_dict.empty()) {
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
