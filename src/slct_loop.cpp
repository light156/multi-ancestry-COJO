#include "include/macojo.h"


// Incremental update for R^{-1} when adding a SNP in the stepwise loop.
int Cohort::calc_R_inv_forward(int append_index)
{   
    if (R_inv_pre.rows() == 0) {
        R_inv_post = MatrixXd::Identity(1, 1);
        R_inv_post_gcta = MatrixXd::Identity(1, 1) / sumstat(append_index, 6);
        return 1;
    }

    RowVectorXd r_temp_row = r.row(append_index);
    RowVectorXd temp_vector = r_temp_row * R_inv_pre;
    if (r_temp_row.dot(temp_vector) > params.collinear) return 0;

    double temp_element = 1.0 / (1.0 - r_temp_row.dot(temp_vector));

    int dim = R_inv_pre.rows();
    R_inv_post = R_inv_pre + temp_element * temp_vector.transpose() * temp_vector;
    R_inv_post.conservativeResize(dim + 1, dim + 1);
    R_inv_post.topRightCorner(dim, 1) = -temp_element * temp_vector.transpose();
    R_inv_post.bottomLeftCorner(1, dim) = -temp_element * temp_vector;
    R_inv_post(dim, dim) = temp_element;

    if (R_inv_post.cwiseAbs().maxCoeff() > params.iter_collinear_threshold) return -1;

    if (params.slct_mode == "GCTA") {
        r_temp_row = r_gcta.row(append_index);
        temp_vector = r_temp_row * R_inv_pre_gcta;
        temp_element = 1.0 / (sumstat(append_index, 6) - r_temp_row.dot(temp_vector));

        R_inv_post_gcta = R_inv_pre_gcta + temp_element * temp_vector.transpose() * temp_vector;
        R_inv_post_gcta.conservativeResize(dim + 1, dim + 1);
        R_inv_post_gcta.topRightCorner(dim, 1) = -temp_element * temp_vector.transpose();
        R_inv_post_gcta.bottomLeftCorner(1, dim) = -temp_element * temp_vector;
        R_inv_post_gcta(dim, dim) = temp_element;
    }

    return 1;
}


// Remove one row/column from R^{-1} when dropping a SNP.
int Cohort::calc_R_inv_backward(int remove_index)
{
    int dim = R_inv_pre.rows();
    if (remove_index < 0 || remove_index >= dim)
        throw invalid_argument("calc_R_inv_backward: remove_index out of range.");

    R_inv_post = R_inv_pre - R_inv_pre.col(remove_index) * R_inv_pre.row(remove_index) / R_inv_pre(remove_index, remove_index);
    remove_row(R_inv_post, remove_index);
    remove_column(R_inv_post, remove_index);

    if (params.slct_mode == "GCTA") {
        R_inv_post_gcta = R_inv_pre_gcta - R_inv_pre_gcta.col(remove_index) * R_inv_pre_gcta.row(remove_index) / R_inv_pre_gcta(remove_index, remove_index);
        remove_row(R_inv_post_gcta, remove_index);
        remove_column(R_inv_post_gcta, remove_index);
    }

    return (R_inv_post.cwiseAbs().maxCoeff() > params.iter_collinear_threshold) ? -1 : 1;
}


// Append correlations between a candidate SNP and active SNPs in-window.
void Cohort::append_r(const vector<char>& active_mask, int append_index, string mode) 
{  
    VectorXd r_temp_vec = VectorXd::Zero(active_mask.size());

    const double left_bp = shared.SNP_pos_ref[append_index] - params.window_size + 1;
    const double right_bp = shared.SNP_pos_ref[append_index] + params.window_size - 1;

    // open interval (left_bp, right_bp)
    auto it_lo = std::lower_bound(shared.bp_order.begin(), shared.bp_order.end(), left_bp,
        [&](int idx, double value) { return shared.SNP_pos_ref[idx] < value; });
    auto it_hi = std::upper_bound(shared.bp_order.begin(), shared.bp_order.end(), right_bp,
        [&](double value, int idx) { return value < shared.SNP_pos_ref[idx]; });

    if (params.if_LD_mode) {
        for (auto it = it_lo; it != it_hi; ++it) {
            int sweep_index = *it;
            if (!active_mask[sweep_index]) continue;
            r_temp_vec(sweep_index) = LD_matrix(sweep_index, append_index);
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(it_hi - it_lo); i++) {
            int sweep_index = *(it_lo + i);
            if (!active_mask[sweep_index]) continue;
            r_temp_vec(sweep_index) = genotype.calc_inner_product(sweep_index, append_index, (mode == "removeNA"));
        }
    }

    append_column(r, r_temp_vec);

    if (mode == "GCTA") {
        double N_new = sumstat(append_index, 4);
        double V_new = sumstat(append_index, 5);
        VectorXd r_temp_vec_gcta = r_temp_vec.array() * sumstat.col(4).min(N_new) * sqrt(sumstat.col(5)) * sqrt(V_new);
        append_column(r_gcta, r_temp_vec_gcta);
    }
}


// Compute conditional effects given candidate SNPs.
void Cohort::calc_cond_effects(const vector<int>& candidate_SNP, string mode)
{   
    if (candidate_SNP.size() == 0) {
        beta = sumstat.col(0);
        beta_var = sumstat.col(1);
        return;
    }
    
    ArrayXXd sumstat_candidate = sumstat(candidate_SNP, all);

    if (mode == "GCTA") {
        VectorXd temp = R_inv_pre_gcta * (sumstat_candidate.col(0) * sumstat_candidate.col(6)).matrix();
        beta = sumstat.col(0) - (r_gcta * temp).array() / sumstat.col(6);
        beta_var = Vp / sumstat.col(6);
    } else {
        VectorXd temp = R_inv_pre * (sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5))).matrix();
        beta = sumstat.col(0) - (r * temp).array() / sqrt(sumstat.col(5));
        beta_var = sumstat.col(1);
    }
}


// Compute joint effects after R^{-1} is ready.
int Cohort::calc_joint_effects(const vector<int>& candidate_SNP, string mode)
{   
    ArrayXXd sumstat_candidate = sumstat(candidate_SNP, all);
    double sigma_J_squared;

    if (candidate_SNP.size() == 1) {
        beta = sumstat_candidate.col(0);
        beta_var = sumstat_candidate.col(1);
        sigma_J_squared = Vp - sumstat_candidate(0, 0) * sumstat_candidate(0, 0) * sumstat_candidate(0, 5);
    } else if (mode == "GCTA") {
        beta = R_inv_post_gcta * (sumstat_candidate.col(0) * sumstat_candidate.col(6)).matrix();
        beta_var = R_inv_post_gcta.diagonal().array() * Vp;
        sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
    } else {
        ArrayXd temp = R_inv_post * (sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5))).matrix();
        beta = temp / sqrt(sumstat_candidate.col(5));

        sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
        beta_var = sigma_J_squared * R_inv_post.diagonal().array() / sumstat_candidate.col(6);    
    }

    double Neff = median(sumstat_candidate.col(4));
    int M = sumstat_candidate.rows();
    R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);

    return (beta_var.minCoeff() > 1e-30) ? 1 : -1;
}


// Inverse-variance meta-analysis across current cohorts.
void MACOJO::inverse_var_meta(ArrayXd &bma, ArrayXd &se2ma, ArrayXd &abs_zma) 
{   
    if (current_list.size() == 1) {
        auto &c = cohorts[current_list[0]];
        bma = c.beta;
        se2ma = c.beta_var;
    } else {
        bma.setZero(cohorts[current_list[0]].beta.rows());
        se2ma.setZero(cohorts[current_list[0]].beta.rows());

        for (int n : current_list) {
            auto &c = cohorts[n]; 
            bma += c.beta / c.beta_var;
            se2ma += 1 / c.beta_var;
        }

        se2ma = 1 / se2ma;
        bma = bma * se2ma;
    }

    abs_zma = abs(bma) / sqrt(se2ma);
}


// Stepwise selection: forward add by conditional p-value, optional backward selection.
void MACOJO::slct_loop() 
{   
    int min_pC_index, max_pJ_index;
    double min_pC, max_pJ;

    string current_SNP_name;
    bool loop_break_indicator = false;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num < params.max_iter_num) {
        auto start = steady_clock::now();

        LOGGER << candidate_SNP.size() << " " << collinear_SNP.size() << " " << backward_SNP.size() << endl;
        
        for (int n : current_list)
            cohorts[n].calc_cond_effects(candidate_SNP, params.slct_mode);

        inverse_var_meta(bC, se2C, abs_zC);
        abs_zC(bad_SNP) = -1;
        abs_zC(candidate_SNP) = -1;
        abs_zC(collinear_SNP) = -1;
        abs_zC(backward_SNP) = -1;
        
        while (true) {
            // forward: res==0: skip, res==1: success, res==-1: numerical issue, remove
            // backward: res==1: success, res!=-1: skip
            int res_flag;

            // select minimal conditional SNP
            min_pC = erfc(abs_zC.maxCoeff(&min_pC_index) / sqrt(2));
            current_SNP_name = shared.goodSNP_table[min_pC_index].first;

            if (min_pC > params.p_value) {
                LOGGER.w("Calculation finished, No more SNPs above p-threshold");
                LOGGER << current_SNP_name << " " << scientific << min_pC << fixed << endl;
                loop_break_indicator = true;
                break;
            }

            // potential new candidate SNP
            candidate_SNP.push_back(min_pC_index);

            for (int n : current_list) {
                res_flag = cohorts[n].calc_R_inv_forward(min_pC_index);
                if (res_flag != 1) break;

                res_flag = cohorts[n].calc_joint_effects(candidate_SNP, params.slct_mode);
                if (res_flag != 1) break;
            }

            // check R2 increment
            for (int n : current_list) {
                if (res_flag == 1 && params.R2_threshold > -1 && cohorts[n].R2 < (1+params.R2_threshold) * cohorts[n].previous_R2) {
                    LOGGER.i("skipped, R2 increment lower than threshold in Cohort " + to_string(n+1), current_SNP_name);
                    res_flag = 0;
                }
            }

            // numerical issue happened during R_inv or joint effect calculation, remove this SNP
            if (res_flag == -1) {
                active_mask[min_pC_index] = 0;
                collinear_SNP.push_back(min_pC_index);
            }

            if (res_flag != 1) {
                candidate_SNP.pop_back();
                abs_zC(min_pC_index) = -1;
                continue;
            }

            // calculate joint p-value and decide whether to accept this SNP
            // must suffice for the first SNP, just avoid complicated logic
            inverse_var_meta(bJ, se2J, abs_zJ);
            max_pJ = erfc(abs_zJ.minCoeff(&max_pJ_index) / sqrt(2));
            // max_pJ = erfc(abs_zJ.tail(candidate_SNP.size()-fixed_candidate_SNP_num).minCoeff(&max_pJ_index) / sqrt(2));
            // max_pJ_index += fixed_candidate_SNP_num;

            // directly accept this SNP and next iteration
            if (max_pJ <= params.p_value) {
                LOGGER << "Conditional min pC: " << scientific << min_pC << fixed << endl;
                LOGGER << "Joint b, se, max pJ: " << bJ(max_pJ_index) << " " 
                    << sqrt(se2J(max_pJ_index)) << " " << scientific << max_pJ << fixed << endl;

                LOGGER.i("New candidate SNP accepted", current_SNP_name);
                active_mask[min_pC_index] = 0;

                for (int n : current_list) {
                    cohorts[n].append_r(active_mask, min_pC_index, params.slct_mode);
                    cohorts[n].save_temp_model();
                }
                break;
            } 
            
            // trivial case for removing current SNP
            if (max_pJ_index < fixed_candidate_SNP_num || max_pJ_index == abs_zJ.rows()-1) {
                // LOGGER.i("removed, because pJ exceeds p-value threshold", current_SNP_name);
                active_mask[min_pC_index] = 0;
                backward_SNP.push_back(min_pC_index);

                candidate_SNP.pop_back();
                abs_zC(min_pC_index) = -1;
                continue;
            } 

            // backward selection must happen now, save current model as backup, then accept the new candidate SNP
            LOGGER << "Conditional min pC: " << scientific << min_pC << fixed << endl;
            LOGGER << "Joint b, se, max pJ: " << bJ(max_pJ_index) << " " 
                << sqrt(se2J(max_pJ_index)) << " " << scientific << max_pJ << fixed << endl;
            LOGGER.i("New candidate SNP added, start backward selection", current_SNP_name);

            candidate_SNP_backup = candidate_SNP;
            backward_SNP_backup = backward_SNP;

            for (int n : current_list) {
                cohorts[n].save_state();
                cohorts[n].append_r(active_mask, min_pC_index, params.slct_mode);
                cohorts[n].save_temp_model();
            }
            
            while (true) {
                LOGGER.i("Candidate SNP removed", shared.goodSNP_table[candidate_SNP[max_pJ_index]].first);

                backward_SNP.push_back(candidate_SNP[max_pJ_index]);
                candidate_SNP.erase(candidate_SNP.begin() + max_pJ_index);

                // don't think this can happen, but just in case
                if (candidate_SNP.size() == 0) {
                    LOGGER.i("Backward selection failed, no candidate SNPs left");
                    res_flag = -1;
                    break;
                }

                for (int n : current_list) {
                    remove_column(cohorts[n].r, max_pJ_index);
                    if (params.slct_mode == "GCTA")
                        remove_column(cohorts[n].r_gcta, max_pJ_index); 

                    // do not need to check collinearity because it is guaranteed in forward step
                    res_flag = cohorts[n].calc_R_inv_backward(max_pJ_index);
                    if (res_flag != 1) {
                        // don't think this can happen, but just in case
                        LOGGER.i("Backward selection failed, collinearity issue after removing " + current_SNP_name);
                        break;
                    }

                    res_flag = cohorts[n].calc_joint_effects(candidate_SNP, params.slct_mode);
                    if (res_flag != 1) {
                        // don't think this can happen, but just in case
                        LOGGER.i("Backward selection failed, joint se too small after removing " + current_SNP_name);
                        break;
                    }

                    cohorts[n].save_temp_model();
                }

                if (res_flag != 1) break;

                inverse_var_meta(bJ, se2J, abs_zJ);
                max_pJ = erfc(abs_zJ.minCoeff(&max_pJ_index) / sqrt(2));
                // max_pJ = erfc(abs_zJ.tail(candidate_SNP.size()-fixed_candidate_SNP_num).minCoeff(&max_pJ_index) / sqrt(2));
                // max_pJ_index += fixed_candidate_SNP_num;
                LOGGER << "Joint b, se, max pJ: " 
                        << bJ(max_pJ_index) << " " << sqrt(se2J(max_pJ_index)) << " " << scientific << max_pJ << fixed << endl;

                if (max_pJ <= params.p_value) break;
                if (max_pJ_index < fixed_candidate_SNP_num) {
                    LOGGER.i("Backward selection failed, fixed candidate SNP exceeds p-value threshold");
                    res_flag = -1;
                    break;
                }
            } 
            
            // check R2 increment after backward selection
            for (int n : current_list) {
                if (res_flag == 1 && params.R2back_threshold > -1 && cohorts[n].R2 < (1+params.R2back_threshold) * cohorts[n].previous_R2) {
                    LOGGER.i("Backward selection failed, adjusted R2 lower than backward threshold in Cohort " + to_string(n+1));
                    res_flag = -1;
                }
            }

            // restore the last successful model before the newest candidate SNP is added
            if (res_flag != 1) {
                LOGGER.i("Restore to the last successful model before adding " + current_SNP_name);
                candidate_SNP = candidate_SNP_backup;
                backward_SNP = backward_SNP_backup;

                for (int n : current_list)
                    cohorts[n].restore_state();

                candidate_SNP.pop_back();
                abs_zC(min_pC_index) = -1;
                continue;
            }

            // successful backward selection, proceed to next iteration
            active_mask[min_pC_index] = 0;
            LOGGER.i("Backward selection successful!");
            break;
        }

        auto end = steady_clock::now();
        LOGGER << "iter " << iter_num++ << " finished, takes " << duration<double>(end-start).count() << " seconds" << endl;
        LOGGER << "--------------------------------" << endl;
    }

    LOGGER.i("associated SNPs selected", candidate_SNP.size());
    LOGGER.i("SNPs eliminated by backward selection", backward_SNP.size());
    LOGGER.i("SNPs filtered by collinearity test", collinear_SNP.size());
}
