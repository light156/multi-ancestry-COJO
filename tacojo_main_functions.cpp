#include "trans_ancestry_cojo.h"


void Cohort::calc_inner_product(const vector<int> &index_list, int single_index, int window_size) 
{   
    r_temp_vec.setZero(index_list.size());

    # pragma omp parallel for 
    for (int i = 0; i < index_list.size(); i++) {
        if (abs(SNP_pos_ref[index_list[i]] - SNP_pos_ref[single_index]) <= window_size)
            r_temp_vec(i) = (X.col(single_index).transpose() * X.col(index_list[i])).value() / (indi_num-1);
    }
}


void Cohort::calc_conditional_effects() 
{   
    MatrixXd temp1, temp2;
    ArrayXd temp3;

    temp1 = sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5));
    temp2.noalias() = R_inv_pre * temp1;
    temp3 = r * temp2;
    conditional_beta = sumstat_screened.col(0) - temp3 / sqrt(sumstat_screened.col(5));
}


bool Cohort::calc_joint_effects(const ArrayXXd &sumstat_temp, bool flag, double iter_colinear_threshold) 
{      
    if (flag && ((abs(R_inv_post.minCoeff()) > iter_colinear_threshold) || (abs(R_inv_post.maxCoeff()) > iter_colinear_threshold)))
        return true;

    VectorXd temp1 = sqrt(sumstat_temp.col(5)) * sumstat_temp.col(0);
    ArrayXd temp2 = R_inv_post * temp1;
    beta = temp2 / sqrt(sumstat_temp.col(5));

    double sigma_J_squared = Vp - (sumstat_temp.col(0) * sumstat_temp.col(5) * beta).sum();
                
    beta_var = sigma_J_squared * R_inv_post.diagonal().array() / sumstat_temp.col(6) ;

    if (flag && beta_var.minCoeff() <= 0)
        return true;

    double Neff = median(sumstat_temp.col(4));
    int M = sumstat_temp.rows();
    R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    return false;
}


void Cohort::calc_R_inv(bool if_fast_inv) {
    if (if_fast_inv) {
        double temp_element = 1 / (1 - r_temp_vec.transpose() * R_inv_pre * r_temp_vec);
        VectorXd temp_vector = R_inv_pre * r_temp_vec;
    
        int dim = R_inv_pre.rows();
        R_inv_post.resize(dim+1, dim+1);
    
        R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
        R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
        R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
        R_inv_post(dim, dim) = temp_element;
    } else {
        int dim = R_pre.rows();
        R_post.resize(dim+1, dim+1);
    
        R_post.block(0, 0, dim, dim) = R_pre;
        R_post.block(0, dim, dim, 1) = r_temp_vec;
        R_post.block(dim, 0, 1, dim) = r_temp_vec.transpose();
        R_post(dim, dim) = 1;

        R_inv_post.noalias() = R_post.llt().solve(MatrixXd::Identity(dim+1, dim+1));
    }
}


void TransAncestryCOJO::inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
    const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge) 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    merge.resize(b_cohort1.size(), 4);

    merge.col(0) = (b_cohort1 * se2_cohort2 + b_cohort2 * se2_cohort1) / (se2_cohort1 + se2_cohort2);
    merge.col(1) = se2_cohort1 * se2_cohort2 / (se2_cohort1 + se2_cohort2);
    merge.col(2) = abs(merge.col(0) / sqrt(merge.col(1)));
    merge.col(3) = erfc(merge.col(2) / sqrt(2));
}


void TransAncestryCOJO::initialize_matrices(Cohort &c) 
{   
    c.calc_inner_product(screened_SNP, max_SNP_index, window_size);

    c.sumstat_candidate = c.sumstat.row(max_SNP_index);
    c.sumstat_screened = c.sumstat;

    c.r = c.r_temp_vec;
    c.R_inv_pre = MatrixXd::Identity(1,1);
    c.R_pre = MatrixXd::Identity(1,1);
    
    c.previous_R2 = 0.0;
    c.output_b = c.sumstat.col(0);
    c.output_se2 = c.sumstat.col(1);
}


void TransAncestryCOJO::initialize_MDISA(Cohort &c) 
{   
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {
        c.calc_inner_product(candidate_SNP, excluded_SNP[i], window_size); 
        
        if (c.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) {
            append_row(c.sumstat_screened, c.sumstat.row(excluded_SNP[i])); 
            append_row(c.r, c.r_temp_vec.transpose());
            
            screened_SNP.push_back(excluded_SNP[i]); 
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }

    for (int i = backward_removed_SNP.size()-1; i >= 0; i--) {
        c.calc_inner_product(candidate_SNP, backward_removed_SNP[i], window_size); 
        
        if (c.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) {
            append_row(c.sumstat_screened, c.sumstat.row(backward_removed_SNP[i])); 
            append_row(c.r, c.r_temp_vec.transpose());

            screened_SNP.push_back(backward_removed_SNP[i]);
        } else {
            excluded_SNP.push_back(backward_removed_SNP[i]); 
        }
    }

    vector<int>().swap(backward_removed_SNP);
    c.R_inv_post = c.R_inv_pre;
}


void TransAncestryCOJO::initialize_backward_selection(Cohort &c, const ArrayXd &pJ)
{   
    c.sumstat_backward = c.sumstat_candidate;

    for (int i = candidate_SNP.size()-1; i >= fixed_candidate_SNP_num; i--) {
        if (pJ(i) > threshold)
            remove_row(c.sumstat_backward, i);
    }

    vector<VectorXd> X_backward;
    vector<int> SNP_pos_backward;
    int total_num = 0;
    
    for (int i = 0; i < candidate_SNP.size(); i++) {
        if (i < fixed_candidate_SNP_num || pJ(i) <= threshold) {
            X_backward.push_back(c.X.col(candidate_SNP[i]));
            SNP_pos_backward.push_back(Cohort::SNP_pos_ref[candidate_SNP[i]]);
            total_num++;
        }
    }

    c.R_post = MatrixXd::Identity(total_num, total_num);

    # pragma omp parallel for 
    for (int i = 1; i < total_num; i++)
        for (int j = 0; j < i; j++)
            if (abs(SNP_pos_backward[i] - SNP_pos_backward[j]) <= window_size) {
                c.R_post(i, j) = (X_backward[i].transpose() * X_backward[j]).value() / (c.indi_num-1);
                c.R_post(j, i) = c.R_post(i, j);
            } 

    vector<VectorXd>().swap(X_backward);
    vector<int>().swap(SNP_pos_backward);
        
    // c.R_inv_post.noalias() = R_post_cohort.inverse();
    c.R_inv_post.noalias() = c.R_post.llt().solve(MatrixXd::Identity(total_num, total_num));
}


void TransAncestryCOJO::remove_new_colinear_SNP(bool cohort1_only, bool cohort2_only) 
{   
    for (int i = screened_SNP.size()-1, last_col = candidate_SNP.size()-1; i >= 0; i--) {
        if ((!cohort2_only && abs(c1.r(i, last_col)) >= colinear_threshold_sqrt) || 
            (!cohort1_only && abs(c2.r(i, last_col)) >= colinear_threshold_sqrt)) {
            if (!cohort2_only) {
                remove_row(c1.sumstat_screened, i);
                remove_row(c1.r, i);
            }

            if (!cohort1_only) {
                remove_row(c2.sumstat_screened, i);
                remove_row(c2.r, i);
            }
            
            if (i != max_SNP_index)
                excluded_SNP.push_back(screened_SNP[i]);
            screened_SNP.erase(screened_SNP.begin()+i);
        }
    }
}   


void TransAncestryCOJO::adjust_SNP_according_to_backward_selection(const ArrayXd &pJ, bool cohort1_only, bool cohort2_only) 
{   
    if (!cohort2_only) {
        remove_row(c1.sumstat_screened, max_SNP_index);
        remove_row(c1.r, max_SNP_index);
    }

    if (!cohort1_only) {
        remove_row(c2.sumstat_screened, max_SNP_index);
        remove_row(c2.r, max_SNP_index);
    }

    screened_SNP.erase(screened_SNP.begin()+max_SNP_index);

    for (int i = candidate_SNP.size()-1; i >= fixed_candidate_SNP_num; i--) {
        if (pJ(i) > threshold) {
            if (!cohort2_only) {
                remove_row(c1.sumstat_candidate, i);
                remove_column(c1.r, i);
            }

            if (!cohort1_only) {
                remove_row(c2.sumstat_candidate, i);
                remove_column(c2.r, i);  
            }

            LOGGER.w(1, "Previous candidate SNP removed", Cohort::commonSNP[candidate_SNP[i]]);
            backward_removed_SNP.push_back(candidate_SNP[i]);
            candidate_SNP.erase(candidate_SNP.begin()+i);
        }
    }
    
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {
        if (!cohort2_only)
            c1.calc_inner_product(candidate_SNP, excluded_SNP[i], window_size);

        if (!cohort1_only)
            c2.calc_inner_product(candidate_SNP, excluded_SNP[i], window_size); 
        
        if ((cohort2_only || c1.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) && 
            (cohort1_only || c2.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt)) {
            if (!cohort2_only) {
                append_row(c1.sumstat_screened, c1.sumstat.row(excluded_SNP[i]));
                append_row(c1.r, c1.r_temp_vec.transpose());
            }

            if (!cohort1_only) {
                append_row(c2.sumstat_screened, c2.sumstat.row(excluded_SNP[i])); 
                append_row(c2.r, c2.r_temp_vec.transpose());
            }
            
            screened_SNP.push_back(excluded_SNP[i]); 
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }
}


void Cohort::save_temp_model(bool if_fast_inv) 
{   
    if (if_fast_inv)
        R_inv_pre = R_inv_post;
    else
        R_pre = R_post;

    previous_R2 = R2;
    output_b = beta;
    output_se2 = beta_var;
}


void TransAncestryCOJO::main_loop(string savename) 
{   
    inverse_var_meta(c1.sumstat.col(0), c2.sumstat.col(0), c1.sumstat.col(1), c2.sumstat.col(1), sumstat_merge);
    
    sumstat_merge.col(2).maxCoeff(&max_SNP_index);
    if (sumstat_merge(max_SNP_index, 3) > threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    LOGGER.i(0, "First SNP", Cohort::commonSNP[max_SNP_index]);
    LOGGER << "--------------------------------" << endl;

    candidate_SNP.push_back(max_SNP_index);
    screened_SNP.resize(Cohort::commonSNP.size());
    iota(screened_SNP.begin(), screened_SNP.end(), 0);

    initialize_matrices(c1);
    initialize_matrices(c2);
    remove_new_colinear_SNP();

    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;
    string max_SNP_name;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " "
            << excluded_SNP.size() << " " << backward_removed_SNP.size() << endl;

        // calculate conditional effects
        c1.calc_conditional_effects();
        c2.calc_conditional_effects();
        inverse_var_meta(c1.conditional_beta, c2.conditional_beta, 
            c1.sumstat_screened.col(1), c2.sumstat_screened.col(1), sumstat_merge);
        
        while (true) {
            // select maximal SNP
            double max_Zabs = sumstat_merge.col(2).maxCoeff(&max_SNP_index);
            if (max_Zabs < 0 || sumstat_merge(max_SNP_index, 3) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }
            
            max_SNP_name = Cohort::commonSNP[screened_SNP[max_SNP_index]];

            // calculate joint effects
            append_row(c1.sumstat_candidate, c1.sumstat_screened.row(max_SNP_index));
            c1.calc_inner_product(candidate_SNP, screened_SNP[max_SNP_index], window_size);
            c1.calc_R_inv(if_fast_inv);
            NA_flag = c1.calc_joint_effects(c1.sumstat_candidate, true, iter_colinear_threshold);

            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            append_row(c2.sumstat_candidate, c2.sumstat_screened.row(max_SNP_index));
            c2.calc_inner_product(candidate_SNP, screened_SNP[max_SNP_index], window_size);
            c2.calc_R_inv(if_fast_inv);
            NA_flag = c2.calc_joint_effects(c2.sumstat_candidate, true, iter_colinear_threshold);

            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                remove_row(c2.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            inverse_var_meta(c1.beta, c2.beta, c1.beta_var, c2.beta_var, sumstat_new_model_joint);

            if ((c1.R2 < (1+R2_incremental_threshold) * c1.previous_R2) || 
                (c2.R2 < (1+R2_incremental_threshold) * c2.previous_R2)) {
                LOGGER.w(1, "R2 increment unsatisfactory", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                remove_row(c2.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            // include new candidate SNP
            candidate_SNP.push_back(screened_SNP[max_SNP_index]);

            c1.calc_inner_product(screened_SNP, screened_SNP[max_SNP_index], window_size);
            append_column(c1.r, c1.r_temp_vec);
            
            c2.calc_inner_product(screened_SNP, screened_SNP[max_SNP_index], window_size);
            append_column(c2.r, c2.r_temp_vec);
            
            if (sumstat_new_model_joint.col(3).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed", max_SNP_name);

                int M = candidate_SNP.size();
                // LOGGER << "Added R vector cohort 1: " << c1.r.row(max_SNP_index) << endl;
                // LOGGER << "Added R vector cohort 2: " << c2.r.row(max_SNP_index) << endl;
                remove_new_colinear_SNP();
                LOGGER << "Added diagonal value cohort 1: " << c1.R_inv_post(M-1, M-1) << endl;
                LOGGER << "Added diagonal value cohort 2: " << c2.R_inv_post(M-1, M-1) << endl;
                LOGGER << "Joint b: " << sumstat_new_model_joint(M-1, 0) << endl;
                LOGGER << "Joint se: " << sqrt(sumstat_new_model_joint(M-1, 1)) << endl;
                LOGGER << "Joint p-value: " << scientific << sumstat_new_model_joint(M-1, 3) << endl;
                LOGGER << "Adjusted R2 for cohort 1: " << fixed << c1.R2 << endl;
                LOGGER << "Adjusted R2 for cohort 2: " << c2.R2 << endl;
                break; 
            }

            // backward selection
            LOGGER << "Backward selection" << endl;

            initialize_backward_selection(c1, sumstat_new_model_joint.col(3));
            initialize_backward_selection(c2, sumstat_new_model_joint.col(3));

            c1.calc_joint_effects(c1.sumstat_backward, false);
            c2.calc_joint_effects(c2.sumstat_backward, false);

            LOGGER << "Adjusted R2 after SNP elimination for cohort 1: " << fixed << c1.R2 << endl;
            LOGGER << "Adjusted R2 after SNP elimination for cohort 2: " << c2.R2 << endl;

            if ((c1.R2 < (1+R2_incremental_threshold_backwards) * c1.previous_R2) || 
                (c2.R2 < (1+R2_incremental_threshold_backwards) * c2.previous_R2)) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", max_SNP_name);
                candidate_SNP.pop_back();
                remove_column(c1.r);
                remove_column(c2.r);
                remove_row(c1.sumstat_candidate);
                remove_row(c2.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", max_SNP_name);
            adjust_SNP_according_to_backward_selection(sumstat_new_model_joint.col(3));
            remove_new_colinear_SNP();
            break;
        }

        // save template model for output
        if (!loop_break_indicator) {
            c1.save_temp_model(if_fast_inv);
            c2.save_temp_model(if_fast_inv);
        }

        LOGGER << "iter " << ++iter_num << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }

    inverse_var_meta(c1.output_b, c2.output_b, c1.output_se2, c2.output_se2, sumstat_merge);
    save_results_main_loop(savename + ".jma.cojo");

    // backup SNPs for Cohort 1
    fixed_candidate_SNP_num = candidate_SNP.size();
    
    vector<int> candidate_SNP_backup = candidate_SNP;
    vector<int> screened_SNP_backup = screened_SNP;
    vector<int> excluded_SNP_backup = excluded_SNP;
    vector<int> backward_removed_SNP_backup = backward_removed_SNP;

    // MDISA Cohort 1
    LOGGER.i(0, "Cohort 1", "MDISA");
    initialize_MDISA(c1);
    MDISA(c1);
    save_results_DISA(c1, savename + ".MDISA.cohort1.jma.cojo");

    // backup SNPs for Cohort 2
    candidate_SNP = candidate_SNP_backup;
    screened_SNP = screened_SNP_backup;
    excluded_SNP = excluded_SNP_backup;
    backward_removed_SNP = backward_removed_SNP_backup;

    // MDISA Cohort 2
    LOGGER.i(0, "Cohort 2", "MDISA");
    initialize_MDISA(c2);
    MDISA(c2);
    save_results_DISA(c2, savename + ".MDISA.cohort2.jma.cojo");

    vector<int>().swap(candidate_SNP_backup);
    vector<int>().swap(screened_SNP_backup);
    vector<int>().swap(excluded_SNP_backup);
    vector<int>().swap(backward_removed_SNP_backup);
}


void TransAncestryCOJO::MDISA(Cohort &c) 
{   
    bool cohort1_only = (addressof(c) == addressof(c1));
    ArrayXd Zabs_temp, p_temp, ZabsJ, pJ;

    if (candidate_SNP.size() == 0) {

        Zabs_temp = abs(c.sumstat.col(0) / c.sumstat.col(1));
        p_temp = erfc(Zabs_temp / sqrt(2));

        Zabs_temp.maxCoeff(&max_SNP_index);
        if (p_temp(max_SNP_index) > threshold)
            LOGGER.e(0, "Input data has no significant SNPs.");

        LOGGER.i(0, "First SNP", Cohort::commonSNP[max_SNP_index]);

        candidate_SNP.push_back(max_SNP_index);
        screened_SNP.resize(Cohort::commonSNP.size());
        iota(screened_SNP.begin(), screened_SNP.end(), 0);

        // impossible to happen, just in case
        vector<int>().swap(excluded_SNP);
        vector<int>().swap(backward_removed_SNP);
        
        initialize_matrices(c);
        remove_new_colinear_SNP(cohort1_only, !cohort1_only);
    } 

    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;
    string max_SNP_name;

    LOGGER << "--------------------------------" << endl;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " "
            << excluded_SNP.size() << " " << backward_removed_SNP.size() << endl;

        // calculate conditional effects
        c.calc_conditional_effects();

        Zabs_temp = abs(c.conditional_beta / sqrt(c.sumstat_screened.col(1)));
        p_temp = erfc(Zabs_temp / sqrt(2));
    
        while (true) {
            // select maximal SNP
            double max_Zabs = Zabs_temp.maxCoeff(&max_SNP_index);
            if (max_Zabs < 0 || p_temp(max_SNP_index) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }
            
            max_SNP_name = Cohort::commonSNP[screened_SNP[max_SNP_index]];

            // calculate joint effects
            append_row(c.sumstat_candidate, c.sumstat_screened.row(max_SNP_index));
            c.calc_inner_product(candidate_SNP, screened_SNP[max_SNP_index], window_size);
            c.calc_R_inv(if_fast_inv);
            NA_flag = c.calc_joint_effects(c.sumstat_candidate, true, iter_colinear_threshold);

            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c.sumstat_candidate);
                Zabs_temp(max_SNP_index) = -1;
                continue;
            }

            if (c.R2 < (1+R2_incremental_threshold) * c.previous_R2) {
                LOGGER.w(1, "R2 increment unsatisfactory", max_SNP_name);
                remove_row(c.sumstat_candidate);
                Zabs_temp(max_SNP_index) = -1;
                continue;
            }

            ZabsJ = abs(c.beta / sqrt(c.beta_var));
            pJ = erfc(ZabsJ / sqrt(2));

            // include new candidate SNP
            candidate_SNP.push_back(screened_SNP[max_SNP_index]);

            c.calc_inner_product(screened_SNP, screened_SNP[max_SNP_index], window_size);
            append_column(c.r, c.r_temp_vec);
            
            if (pJ.bottomRows(candidate_SNP.size()-fixed_candidate_SNP_num).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed", max_SNP_name);

                int M = candidate_SNP.size();
                // LOGGER << "Added R vector: " << c.r.row(max_SNP_index) << endl;
                remove_new_colinear_SNP(cohort1_only, !cohort1_only);
                LOGGER << "Added diagonal value: " << c.R_inv_post(M-1, M-1) << endl;
                LOGGER << "Joint b: " << c.beta(M-1) << endl;
                LOGGER << "Joint se: " << sqrt(c.beta_var(M-1)) << endl;
                LOGGER << "Joint p-value: " << scientific << pJ(M-1) << endl;
                LOGGER << "Adjusted R2: " << fixed << c.R2 << endl;
                break; 
            }

            // backward selection
            LOGGER << "Backward selection" << endl;

            initialize_backward_selection(c, pJ);
            c.calc_joint_effects(c.sumstat_backward, false);

            LOGGER << "Adjusted R2 after SNP elimination: " << fixed << c.R2 << endl;

            if (c.R2 < (1+R2_incremental_threshold_backwards) * c.previous_R2) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", max_SNP_name);
                candidate_SNP.pop_back();
                remove_column(c.r);
                remove_row(c.sumstat_candidate);
                Zabs_temp(max_SNP_index) = -1;
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", max_SNP_name);
            adjust_SNP_according_to_backward_selection(pJ, cohort1_only, !cohort1_only);
            remove_new_colinear_SNP(cohort1_only, !cohort1_only);
            break;
        }

        // save template model for output
        if (!loop_break_indicator)
            c.save_temp_model(if_fast_inv);

        LOGGER << "iter " << ++iter_num << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }
}


void TransAncestryCOJO::save_results_main_loop(string filepath) 
{   
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    string headers[29] = {"SNP", "A1", "A2", 
        "freq.x", "b.x", "se.x", "p.x", "N.x", "D.x", "V.x", "bJ.x", "seJ.x", "zJ.x", "pJ.x", 
        "freq.y", "b.y", "se.y", "p.y", "N.y", "D.y", "V.y", "bJ.y", "seJ.y", "zJ.y", "pJ.y", 
        "seJ.ma", "bJ.ma", "zJ", "pJ"};
    
    for (int i = 0; i < 28; i++) 
        jmaCOJO << headers[i] << "\t";
        
    jmaCOJO << headers[28] << "\n";

    int index = 0;
    map<string, int> SNP_order_pair;
 
    // Inserting element in pair vector to keep track of previous indexes
    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, index++) 
        SNP_order_pair.insert(make_pair(Cohort::commonSNP[*iter], index));
    
    ArrayXd temp_row;
    ArrayXd bJx = c1.output_b, bJy = c2.output_b;
    ArrayXd seJx = sqrt(c1.output_se2), seJy = sqrt(c2.output_se2);
    ArrayXd zJx = bJx / seJx, zJy = bJy / seJy;
    ArrayXd pJx = erfc(abs(zJx)/sqrt(2)), pJy = erfc(abs(zJy)/sqrt(2));
    
    ArrayXd bJma = sumstat_merge.col(0);
    ArrayXd seJma = sqrt(sumstat_merge.col(1));
    ArrayXd zJma = bJma / seJma;
    ArrayXd pJma = erfc(abs(zJma)/sqrt(2));
    
    jmaCOJO.precision(12);

    for (auto iter = SNP_order_pair.begin(); iter != SNP_order_pair.end(); iter++) {
        index = iter->second;
        jmaCOJO << iter->first << "\t" << Cohort::A1_ref[candidate_SNP[index]] << "\t" << Cohort::A2_ref[candidate_SNP[index]] << "\t";

        // cols 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
        temp_row = c1.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t" << temp_row(6) << "\t" << temp_row(5) << "\t";
        jmaCOJO << bJx(index) << "\t" << seJx(index) << "\t" << zJx(index) << "\t" << pJx(index) << "\t";
            
        temp_row = c2.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t" << temp_row(6) << "\t" << temp_row(5) << "\t";
        jmaCOJO << bJy(index) << "\t" << seJy(index) << "\t" << zJy(index) << "\t" << pJy(index) << "\t";
        
        jmaCOJO << bJma(index) << "\t" << seJma(index) << "\t" << zJma(index) << "\t" << pJma(index) << "\n";
    }
    jmaCOJO.clear();
    jmaCOJO.close();

    map<string, int>().swap(SNP_order_pair);
    LOGGER.i(0, "Results saved into [" + filepath + "]\n", "Finished");
}


void TransAncestryCOJO::save_results_DISA(Cohort &c, string filepath) 
{   
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    string headers[14] = {"SNP", "A1", "A2", "freq", "b", "se", "p", "N", "D", "V", "bJ", "seJ", "zJ", "pJ"};
        
    for (int i = 0; i < 13; i++) 
        jmaCOJO << headers[i] << "\t";
        
    jmaCOJO << headers[13] << "\n";

    int index = 0;
    map<string, int> SNP_order_pair;
 
    // Inserting element in pair vector to keep track of previous indexes
    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, index++) 
        SNP_order_pair.insert(make_pair(Cohort::commonSNP[*iter], index));
    
    ArrayXd temp_row;
    ArrayXd bJ = c.output_b;
    ArrayXd seJ = sqrt(c.output_se2);
    ArrayXd zJ = bJ / seJ;
    ArrayXd pJ = erfc(abs(zJ)/sqrt(2));

    jmaCOJO.precision(12);

    for (auto iter = SNP_order_pair.begin(); iter != SNP_order_pair.end(); iter++) {
        index = iter->second;
        jmaCOJO << iter->first << "\t" << Cohort::A1_ref[candidate_SNP[index]] << "\t" << Cohort::A2_ref[candidate_SNP[index]] << "\t";

        // cols 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
        temp_row = c.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t" << temp_row(6) << "\t" << temp_row(5) << "\t";

        jmaCOJO << bJ(index) << "\t" << seJ(index) << "\t" << zJ(index) << "\t" << pJ(index) << "\n";
    }
    jmaCOJO.clear();
    jmaCOJO.close();

    map<string, int>().swap(SNP_order_pair);
    LOGGER.i(0, "Results saved into [" + filepath + "]\n", "Finished");
}


bool TransAncestryCOJO::initialize_hyperparameters(int argc, char** argv) 
{
    int temp_num = 0;

    while (temp_num < argc) {
        if (strcmp(argv[temp_num], "-colinear") == 0 && temp_num+1 < argc) {
            colinear_threshold = atof(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-R2") == 0 && temp_num+1 < argc) {
            R2_incremental_threshold = atof(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-R2back") == 0 && temp_num+1 < argc) {
            R2_incremental_threshold_backwards = atof(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-iter_num") == 0 && temp_num+1 < argc) {
            max_iter_num = atoi(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-window") == 0 && temp_num+1 < argc && atoi(argv[temp_num+1]) == -1) {
            window_size = INT_MAX;
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-window") == 0 && temp_num+1 < argc && atoi(argv[temp_num+1]) > 0) {
            window_size = atoi(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "--no_fast_inv") == 0) {
            if_fast_inv = false;
            temp_num += 1;
        } else return false;
    }

    colinear_threshold_sqrt = sqrt(colinear_threshold);
    iter_colinear_threshold = 1 / (1 - colinear_threshold);

    LOGGER << "Threshold: 5e-8" << endl;
    LOGGER << "Colinearity threshold: " << colinear_threshold << endl;
    LOGGER << "R2 incremental threshold: " << R2_incremental_threshold << endl;
    LOGGER << "R2 incremental threshold backwards: " << R2_incremental_threshold_backwards << endl;
    LOGGER << "Maximal iteration number: " << max_iter_num << endl;
    LOGGER << "SNP position window (+/-): " << window_size << endl;
    if (if_fast_inv)
        LOGGER << "Use fast matrix inversion" << endl;
    else
        LOGGER << "Use normal matrix inversion" << endl;

    LOGGER << "--------------------------------" << endl;
    return true;
}