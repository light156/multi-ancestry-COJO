#include "macojo.h"


bool Cohort::calc_R_inv_forward(int append_index)
{   
    if (R_inv_pre.rows() == 0) {
        R_inv_post = MatrixXd::Identity(1, 1);
        R_inv_post_gcta = MatrixXd::Identity(1, 1) / sumstat(append_index, 6);
        return true;
    }

    RowVectorXd r_temp_row = r.row(append_index);
    RowVectorXd temp_vector = r_temp_row * R_inv_pre;
    if (r_temp_row.dot(temp_vector) > params.collinear) return false;

    double temp_element = 1.0 / (1.0 - r_temp_row.dot(temp_vector));

    int dim = R_inv_pre.rows();
    R_inv_post = R_inv_pre + temp_element * temp_vector.transpose() * temp_vector;
    R_inv_post.conservativeResize(dim + 1, dim + 1);
    R_inv_post.topRightCorner(dim, 1) = -temp_element * temp_vector.transpose();
    R_inv_post.bottomLeftCorner(1, dim) = -temp_element * temp_vector;
    R_inv_post(dim, dim) = temp_element;

    if (R_inv_post.cwiseAbs().maxCoeff() > params.iter_collinear_threshold)
        return false;

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

    return true;
}


// remove one row and one column
// no need to check colinearity becuase it is guaranteed in forward step
void Cohort::calc_R_inv_backward(int remove_index)
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
}


void Cohort::append_r(const vector<int>& SNP_list, int append_index, string mode) 
{  
    VectorXd r_temp_vec = VectorXd::Zero(shared.total_SNP_num);
    r_temp_vec.setZero(shared.total_SNP_num);

    if (params.if_LD_mode) {
        for (size_t i = 0; i < SNP_list.size(); i++) {
            int sweep_index = SNP_list[i];
            if (shared.chr_ref[sweep_index] == shared.chr_ref[append_index] &&
                abs(shared.SNP_pos_ref[sweep_index] - shared.SNP_pos_ref[append_index]) < params.window_size)
                r_temp_vec(sweep_index) = LD_matrix(sweep_index, append_index);
        }
    } else {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < SNP_list.size(); i++) {
            int sweep_index = SNP_list[i];
            if (shared.chr_ref[sweep_index] == shared.chr_ref[append_index] &&
                abs(shared.SNP_pos_ref[sweep_index] - shared.SNP_pos_ref[append_index]) < params.window_size) {
                r_temp_vec(sweep_index) = genotype.calc_inner_product(sweep_index, append_index, (mode == "removeNA"));
            }  
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


bool Cohort::calc_R_inv_from_SNP_list(const vector<int> &SNP_list, string mode) 
{   
    int total_num = SNP_list.size();
    R_post = MatrixXd::Identity(total_num, total_num);

    for (int i = 0; i < total_num; i++) {
        for (int j = 0; j < i; j++) {
            if (shared.chr_ref[SNP_list[i]] == shared.chr_ref[SNP_list[j]] &&
                abs(shared.SNP_pos_ref[SNP_list[i]] - shared.SNP_pos_ref[SNP_list[j]]) < params.window_size) {
                
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

    return R_inv_post.cwiseAbs().maxCoeff() < params.iter_collinear_threshold;
}


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
    }
    else {
        VectorXd temp = R_inv_pre * (sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5))).matrix();
        beta = sumstat.col(0) - (r * temp).array() / sqrt(sumstat.col(5));
        beta_var = sumstat.col(1);
    }
}


bool Cohort::calc_joint_effects(const vector<int>& candidate_SNP, string mode)
{   
    ArrayXXd sumstat_candidate = sumstat(candidate_SNP, all);
    
    if (mode == "GCTA") {
        if (candidate_SNP.size() == 1) {
            beta = sumstat_candidate.col(0);
            beta_var = sumstat_candidate.col(1);
        } else {
            beta = R_inv_post_gcta * (sumstat_candidate.col(0) * sumstat_candidate.col(6)).matrix();
            beta_var = R_inv_post_gcta.diagonal().array() * Vp;
        }
    } else {
        double sigma_J_squared;

        if (candidate_SNP.size() == 1) {
            beta = sumstat_candidate.col(0);
            beta_var = sumstat_candidate.col(1);
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
    }

    return beta_var.minCoeff() > 1e-30;
}


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


void MACOJO::slct_loop() 
{   
    int min_pC_index, max_pJ_index;
    double min_pC, max_pJ;

    string current_SNP_name;
    bool loop_break_indicator = false, success_flag, backward_success_flag;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num < params.max_iter_num && screened_SNP.size() > 0) {
        auto start = chrono::steady_clock::now();

        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " " << collinear_SNP.size() << " " << backward_SNP.size() << endl;
        
        for (int n : current_list)
            cohorts[n].calc_cond_effects(candidate_SNP, params.slct_mode);

        inverse_var_meta(bC, se2C, abs_zC);
        abs_zC(candidate_SNP) = -1;
        abs_zC(collinear_SNP) = -1;
        abs_zC(backward_SNP) = -1;
        abs_zC(bad_SNP) = -1;
        
        while (true) {
            success_flag = true;
            backward_success_flag = true;

            // select minimal conditional SNP
            min_pC = erfc(abs_zC.maxCoeff(&min_pC_index) / sqrt(2));
            current_SNP_name = shared.SNP_ref[min_pC_index];

            if (min_pC > params.p) {
                LOGGER.w(0, "Calculation finished, No more SNPs above p-threshold");
                LOGGER << current_SNP_name << " " << scientific << min_pC << fixed << endl;
                loop_break_indicator = true;
                break;
            }

            // potential new candidate SNP
            candidate_SNP.push_back(min_pC_index);

            // calculate R_inv and joint effects for all cohorts
            for (int n : current_list) {
                if (!cohorts[n].calc_R_inv_forward(min_pC_index)) {
                    // LOGGER.w(1, "skipped, colinearity exceeds threshold in Cohort " + to_string(n+1), current_SNP_name);
                    success_flag = false;
                    break;
                }

                if (!cohorts[n].calc_joint_effects(candidate_SNP, params.slct_mode)) {
                    LOGGER.w(1, "skipped, joint se too small in Cohort " + to_string(n+1), current_SNP_name);
                    success_flag = false;
                    break;
                }
            }

            if (!success_flag) {
                screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));
                collinear_SNP.push_back(min_pC_index);

                candidate_SNP.pop_back();
                abs_zC(min_pC_index) = -1;
                continue;
            }

            // check R2 increment, which does not need to be done in gcta-COJO
            if (params.slct_mode != "GCTA") {
                for (int n : current_list) {
                    if (cohorts[n].R2 < (1+params.R2_threshold) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "skipped, R2 increment lower than threshold in Cohort " + to_string(n+1), current_SNP_name);
                        success_flag = false;
                    }
                }

                if (!success_flag) {
                    candidate_SNP.pop_back();
                    abs_zC(min_pC_index) = -1;
                    continue;
                }
            }
            
            LOGGER << "Conditional min pC: " << scientific << min_pC << fixed << endl;

            // calculate joint p-value and decide whether to accept this SNP
            // must suffice for the first SNP, just avoid complicated logic
            inverse_var_meta(bJ, se2J, abs_zJ);
            max_pJ = erfc(abs_zJ.minCoeff(&max_pJ_index) / sqrt(2));
            LOGGER << "Joint b, se, max pJ: " 
                    << bJ(max_pJ_index) << " " << sqrt(se2J(max_pJ_index)) << " " << scientific << max_pJ << fixed << endl;

            // directly accept this SNP and next iteration
            if (max_pJ <= params.p) {
                LOGGER.i(0, "New candidate SNP accepted", current_SNP_name);
                screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));

                for (int n : current_list) {
                    cohorts[n].append_r(screened_SNP, min_pC_index, params.slct_mode);
                    cohorts[n].save_temp_model();
                }
                break;
            } 
            
            // trivial case for removing current SNP
            if (max_pJ_index < fixed_candidate_SNP_num || max_pJ_index == abs_zJ.rows()-1) {
                LOGGER.i(0, "removed, because pJ exceeds p-value threshold", current_SNP_name);
                screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));
                backward_SNP.push_back(min_pC_index);

                candidate_SNP.pop_back();
                abs_zC(min_pC_index) = -1;
                continue;
            } 

            // backward selection must happen now, save current model as backup, then accept the new candidate SNP
            LOGGER.i(0, "New candidate SNP added, start backward selection", current_SNP_name);

            candidate_SNP_backup = candidate_SNP;
            backward_SNP_backup = backward_SNP;

            for (int n : current_list) {
                cohorts[n].save_state();
                cohorts[n].append_r(screened_SNP, min_pC_index, params.slct_mode);
                cohorts[n].save_temp_model();
            }

            // backward selection
            while (true) {
                LOGGER.i(0, "Candidate SNP removed", shared.SNP_ref[candidate_SNP[max_pJ_index]]);

                backward_SNP.push_back(candidate_SNP[max_pJ_index]);
                candidate_SNP.erase(candidate_SNP.begin() + max_pJ_index);

                // don't think this can happen, but just in case
                if (candidate_SNP.size() == 0) {
                    LOGGER.i(0, "Backward selection failed, no candidate SNPs left");
                    backward_success_flag = false;
                    break;
                }

                for (int n : current_list) {
                    remove_column(cohorts[n].r, max_pJ_index);
                    if (params.slct_mode == "GCTA")
                        remove_column(cohorts[n].r_gcta, max_pJ_index); 

                    // do not need to check colinearity because it is guaranteed in forward step
                    cohorts[n].calc_R_inv_backward(max_pJ_index);

                    if (!cohorts[n].calc_joint_effects(candidate_SNP, params.slct_mode)) {
                        // don't think this can happen, but just in case
                        LOGGER.i(0, "Backward selection failed, joint se too small after removing " + current_SNP_name);
                        backward_success_flag = false;
                        break;
                    }

                    cohorts[n].save_temp_model();
                }

                if (!backward_success_flag) break;

                inverse_var_meta(bJ, se2J, abs_zJ);
                max_pJ = erfc(abs_zJ.minCoeff(&max_pJ_index) / sqrt(2));
                LOGGER << "Joint b, se, max pJ: " 
                        << bJ(max_pJ_index) << " " << sqrt(se2J(max_pJ_index)) << " " << scientific << max_pJ << fixed << endl;

                if (max_pJ <= params.p) break;
                if (max_pJ_index < fixed_candidate_SNP_num) {
                    LOGGER.i(0, "Backward selection failed, fixed candidate SNP exceeds p-value threshold");
                    backward_success_flag = false;
                    break;
                }
            } 
            
            // check R2 increment after backward selection, which does not need to be done in gcta-COJO
            if (backward_success_flag && params.slct_mode != "GCTA") {
                for (int n : current_list) {
                    if (cohorts[n].R2 < (1+params.R2back_threshold) * cohorts[n].previous_R2) {
                        LOGGER.w(1, "Backward selection failed, adjusted R2 lower than backward threshold", current_SNP_name);
                        backward_success_flag = false;
                    }
                }
            }

            // restore the last successful model before the newest candidate SNP is added
            if (!backward_success_flag) {
                LOGGER.i(0, "Restore to the last successful model before adding " + current_SNP_name);
                candidate_SNP = candidate_SNP_backup;
                backward_SNP = backward_SNP_backup;

                for (int n : current_list)
                    cohorts[n].restore_state();

                candidate_SNP.pop_back();
                abs_zC(min_pC_index) = -1;
                continue;
            }

            // successful backward selection, remove index from screened_SNP and proceed to next iteration
            screened_SNP.erase(find(screened_SNP.begin(), screened_SNP.end(), min_pC_index));
            LOGGER.i(0, "Backward selection successful!");
            break;
        }

        auto end = chrono::steady_clock::now();
        LOGGER << "iter " << iter_num++ << " finished, takes " 
                << chrono::duration<double>(end-start).count() << " seconds" << endl;
        LOGGER << "--------------------------------" << endl;
    }

    LOGGER << candidate_SNP.size() << " associated SNPs selected, and " 
            << backward_SNP.size() << " SNPs eliminated by backward selection" << endl << endl;
}
