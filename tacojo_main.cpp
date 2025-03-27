#include "trans_ancestry_cojo.h"


void TransAncestryCOJO::read_files(string cojoFile1, string cojoFile2, string PLINK1, string PLINK2) 
{   
    clock_t tStart;
    
    tStart = clock();

    // read cojofiles and get rough common SNPs    
    c1.read_sumstat(cojoFile1);
    c2.read_sumstat(cojoFile2);
    
    set_intersection(c1.SNP.begin(), c1.SNP.end(), 
        c2.SNP.begin(), c2.SNP.end(), back_inserter(Cohort::commonSNP));
    
    // exclude rare SNP
    int temp_index = 0;
    for (auto iter = Cohort::commonSNP.begin(); iter != Cohort::commonSNP.end(); iter++) {
        if (abs(c1.freq[c1.SNP_index[*iter]]-0.5) <= 0.49 || abs(c2.freq[c2.SNP_index[*iter]]-0.5) <= 0.49) {
            Cohort::sumstat_commonSNP_index.insert(Cohort::sumstat_commonSNP_index.end(), pair<string, int> (*iter, temp_index));
            temp_index++;
        }
    }

    vector<string>().swap(Cohort::commonSNP);

    // set bim file 1 as reference, compare and get common SNPs across 4 files
    Cohort::A1_ref.resize(temp_index);
    Cohort::A2_ref.resize(temp_index);
    Cohort::SNP_pos_ref.resize(temp_index);

    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

    c1.read_PLINK(PLINK1, true);
    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();
    c2.read_PLINK(PLINK2, false);
    
    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

    set_intersection(c1.included_SNP_PLINK.begin(), c1.included_SNP_PLINK.end(), 
        c2.included_SNP_PLINK.begin(), c2.included_SNP_PLINK.end(), back_inserter(Cohort::commonSNP));
    
    if (Cohort::commonSNP.size() == 0)
        LOGGER.e(0, "Input data has no common SNPs.");

    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();

    LOGGER.i(0, "common SNPs included for two cohorts", to_string(Cohort::commonSNP.size()));
    LOGGER.i(0, "individuals in Cohort 1", to_string(c1.indi_num));
    LOGGER.i(0, "individuals in Cohort 2", to_string(c2.indi_num));

    // initialize sumstat and X matrices, calculate Vp
    c1.generate_sumstat_and_X();
    c2.generate_sumstat_and_X();
    LOGGER << "Vp: " << c1.Vp << " " << c2.Vp << endl << endl;
    
    // clean ref A1, A2 and SNP position lists
    vector<string> A1_ref_clean, A2_ref_clean;
    vector<int> SNP_pos_ref_clean;

    for (auto iter = Cohort::commonSNP.begin(); iter != Cohort::commonSNP.end(); iter++) {
        temp_index = Cohort::sumstat_commonSNP_index[*iter];
        A1_ref_clean.push_back(Cohort::A1_ref[temp_index]);
        A2_ref_clean.push_back(Cohort::A2_ref[temp_index]);
        SNP_pos_ref_clean.push_back(Cohort::SNP_pos_ref[temp_index]);
    }
    
    A1_ref_clean.swap(Cohort::A1_ref);
    A2_ref_clean.swap(Cohort::A2_ref);
    SNP_pos_ref_clean.swap(Cohort::SNP_pos_ref);

    vector<string>().swap(A1_ref_clean);
    vector<string>().swap(A2_ref_clean);
    vector<int>().swap(SNP_pos_ref_clean);
    map<string, int>().swap(Cohort::sumstat_commonSNP_index);

    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}


void TransAncestryCOJO::initialize_hyperparameters()
{   
    threshold = 5e-8;
    window_size = 1500000;

    colinear_threshold = 0.8;
    colinear_threshold_sqrt = sqrt(colinear_threshold);
    iter_colinear_threshold = 1 / (1 - colinear_threshold);
    
    R2_incremental_threshold = 0.0;
    R2_incremental_threshold_backwards = -0.5;
    max_iter_num = 1000;  
}
    

void TransAncestryCOJO::inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
    const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge) 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    merge.resize(b_cohort1.size(), 4);

    merge.col(0) = (b_cohort1 * se2_cohort2 + b_cohort2 * se2_cohort1) / (se2_cohort1 + se2_cohort2);
    merge.col(1) = se2_cohort1 * se2_cohort2 / (se2_cohort1 + se2_cohort2);
    merge.col(2) = abs(merge.col(0) / sqrt(merge.col(1)));
    merge.col(3) = erfc(merge.col(2)/sqrt(2));
}


void TransAncestryCOJO::initialize_main_loop(Cohort &c) 
{   
    c.calc_inner_product(screened_SNP, max_SNP_index, window_size);

    c.sumstat_candidate = c.sumstat.row(max_SNP_index);
    c.sumstat_screened = c.sumstat;
    c.r = c.r_temp_vec;
    c.R_inv_pre = MatrixXd::Identity(1,1);
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


void TransAncestryCOJO::initialize_backward_selection(Cohort &c)
{   
    c.sumstat_backward = c.sumstat_candidate;

    for (int i = candidate_SNP.size()-1; i >= fixed_candidate_SNP_num; i--) {
        if (sumstat_new_model_joint(i, 3) > threshold)
            remove_row(c.sumstat_backward, i);
    }

    vector<VectorXd> X_backward;
    vector<int> SNP_pos_backward;
    int total_num = fixed_candidate_SNP_num;
    
    for (int i = 0; i < candidate_SNP.size(); i++) {
        if (i < fixed_candidate_SNP_num || sumstat_new_model_joint(i, 3) <= threshold) {
            X_backward.push_back(c.X.col(candidate_SNP[i]));
            SNP_pos_backward.push_back(Cohort::SNP_pos_ref[candidate_SNP[i]]);
            total_num++;
        }
    }

    MatrixXd R_post_cohort = MatrixXd::Identity(total_num, total_num);

    # pragma omp parallel for 
    for (int i = 1; i < total_num; i++)
        for (int j = 0; j < i; j++)
            if (abs(SNP_pos_backward[i] - SNP_pos_backward[j]) <= window_size) {
                R_post_cohort(i, j) = (X_backward[i].transpose() * X_backward[j]).value() / (c.indi_num-1);
                R_post_cohort(j, i) = R_post_cohort(i, j);
            } 

    vector<VectorXd>().swap(X_backward);
    vector<int>().swap(SNP_pos_backward);
        
    c.R_inv_post.noalias() = R_post_cohort.inverse();
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


void TransAncestryCOJO::adjust_SNP_according_to_backward_selection(bool cohort1_only, bool cohort2_only) 
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
        if (sumstat_new_model_joint(i, 3) > threshold) {
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


void TransAncestryCOJO::main_loop(string savename) 
{   
    initialize_hyperparameters();

    output_b_cohort1 = c1.sumstat.col(0);
    output_b_cohort2 = c2.sumstat.col(0);
    output_se2_cohort1 = c1.sumstat.col(1);
    output_se2_cohort2 = c2.sumstat.col(1);

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);

    sumstat_merge.col(2).maxCoeff(&max_SNP_index);
    if (sumstat_merge(max_SNP_index, 3) > threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    LOGGER.i(0, "First SNP", Cohort::commonSNP[max_SNP_index]);
    LOGGER << "--------------------------------" << endl;

    candidate_SNP.push_back(max_SNP_index);
    screened_SNP.resize(Cohort::commonSNP.size());
    iota(screened_SNP.begin(), screened_SNP.end(), 0);

    initialize_main_loop(c1);
    initialize_main_loop(c2);
    remove_new_colinear_SNP();

    double previous_R2_cohort1 = 0.0, previous_R2_cohort2 = 0.0;
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
            c1.calc_R_inv_fast();
            NA_flag = c1.calc_joint_effects(c1.sumstat_candidate, true, iter_colinear_threshold);

            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            append_row(c2.sumstat_candidate, c2.sumstat_screened.row(max_SNP_index));
            c2.calc_inner_product(candidate_SNP, screened_SNP[max_SNP_index], window_size);
            c2.calc_R_inv_fast();
            NA_flag = c2.calc_joint_effects(c2.sumstat_candidate, true, iter_colinear_threshold);

            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c1.sumstat_candidate);
                remove_row(c2.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }

            inverse_var_meta(c1.beta, c2.beta, c1.beta_var, c2.beta_var, sumstat_new_model_joint);

            if ((c1.R2 < (1+R2_incremental_threshold) * previous_R2_cohort1) || 
                (c2.R2 < (1+R2_incremental_threshold) * previous_R2_cohort2)) {
                // LOGGER.w(1, "R2 increment unsatisfactory", max_SNP_name);
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

            initialize_backward_selection(c1);
            initialize_backward_selection(c2);

            c1.calc_joint_effects(c1.sumstat_backward, false);
            c2.calc_joint_effects(c2.sumstat_backward, false);

            LOGGER << "Adjusted R2 after SNP elimination for cohort 1: " << fixed << c1.R2 << endl;
            LOGGER << "Adjusted R2 after SNP elimination for cohort 2: " << c2.R2 << endl;

            if ((c1.R2 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort1) || 
                (c2.R2 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort2)) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", max_SNP_name);
                candidate_SNP.pop_back();
                remove_column(c1.r);
                remove_row(c1.sumstat_candidate);
                remove_column(c2.r);
                remove_row(c2.sumstat_candidate);
                sumstat_merge(max_SNP_index, 2) = -1;
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", max_SNP_name);
            adjust_SNP_according_to_backward_selection();
            remove_new_colinear_SNP();
            break;
        }

        // save template model for output
        if (!loop_break_indicator) {
            previous_R2_cohort1 = c1.R2;
            previous_R2_cohort2 = c2.R2;

            c1.R_inv_pre = c1.R_inv_post;
            c2.R_inv_pre = c2.R_inv_post;

            output_b_cohort1 = c1.beta;
            output_b_cohort2 = c2.beta;
            output_se2_cohort1 = c1.beta_var;
            output_se2_cohort2 = c2.beta_var; 
        }

        LOGGER << "iter " << ++iter_num << " finished" << endl;
        LOGGER << "--------------------------------" << endl;
    }

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);
    save_results_main_loop(savename + ".jma.cojo");

    // backup SNPs for Cohort 1
    fixed_candidate_SNP_num = candidate_SNP.size();
    
    vector<int> candidate_SNP_backup = candidate_SNP;
    vector<int> screened_SNP_backup = screened_SNP;
    vector<int> excluded_SNP_backup = excluded_SNP;
    vector<int> backward_removed_SNP_backup = backward_removed_SNP;

    // MDISA Cohort 1
    LOGGER.i(0, "MDISA", "Cohort 1");
    initialize_MDISA(c1);
    MDISA(c1);
    save_results_DISA(c1, savename + ".MDISA.cohort1.jma.cojo");

    // backup SNPs for Cohort 2
    candidate_SNP = candidate_SNP_backup;
    screened_SNP = screened_SNP_backup;
    excluded_SNP = excluded_SNP_backup;
    backward_removed_SNP = backward_removed_SNP_backup;

    // MDISA Cohort 2
    LOGGER.i(0, "MDISA", "Cohort 2");
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

    c.calc_joint_effects(c.sumstat_candidate, false);
    c.output_b = c.beta;
    c.output_se2 = c.beta_var;

    double previous_R2 = c.R2;
    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;
    string max_SNP_name;

    ArrayXd Z_temp, p_temp, ZJ, pJ;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        LOGGER << candidate_SNP.size() << " " << screened_SNP.size() << " "
            << excluded_SNP.size() << " " << backward_removed_SNP.size() << endl;

        // calculate conditional effects
        c.calc_conditional_effects();

        Z_temp = abs(c.conditional_beta / sqrt(c.sumstat_screened.col(1)));
        p_temp = erfc(Z_temp/sqrt(2));
    
        while (true) {
            // select maximal SNP
            double max_Zabs = Z_temp.maxCoeff(&max_SNP_index);
            if (max_Zabs < 0 || p_temp(max_SNP_index) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }
            
            max_SNP_name = Cohort::commonSNP[screened_SNP[max_SNP_index]];\

            // calculate joint effects
            append_row(c.sumstat_candidate, c.sumstat_screened.row(max_SNP_index));
            c.calc_inner_product(candidate_SNP, screened_SNP[max_SNP_index], window_size);
            c.calc_R_inv_fast();
            NA_flag = c.calc_joint_effects(c.sumstat_candidate, true, iter_colinear_threshold);

            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", max_SNP_name);
                remove_row(c.sumstat_candidate);
                Z_temp(max_SNP_index) = -1;
                continue;
            }

            if (c.R2 < (1+R2_incremental_threshold) * previous_R2) {
                // LOGGER.w(1, "R2 increment unsatisfactory", max_SNP_name);
                remove_row(c.sumstat_candidate);
                Z_temp(max_SNP_index) = -1;
                continue;
            }

            ZJ = abs(c.beta / sqrt(c.beta_var));
            pJ = erfc(ZJ/sqrt(2));

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

            initialize_backward_selection(c);
            c.calc_joint_effects(c.sumstat_backward, false);

            LOGGER << "Adjusted R2 after SNP elimination: " << fixed << c.R2 << endl;

            if (c.R2 < (1+R2_incremental_threshold_backwards)*previous_R2) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", max_SNP_name);
                candidate_SNP.pop_back();
                remove_column(c.r);
                remove_row(c.sumstat_candidate);
                Z_temp(max_SNP_index) = -1;
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", max_SNP_name);
            adjust_SNP_according_to_backward_selection(cohort1_only, !cohort1_only);
            remove_new_colinear_SNP(cohort1_only, !cohort1_only);
            break;
        }

        // save template model for output
        if (!loop_break_indicator) {
            previous_R2 = c.R2;
            c.R_inv_pre = c.R_inv_post;

            c.output_b = c.beta;
            c.output_se2 = c.beta_var;
        }

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
    ArrayXd bJx = output_b_cohort1, bJy = output_b_cohort2;
    ArrayXd seJx = sqrt(output_se2_cohort1), seJy = sqrt(output_se2_cohort2);
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