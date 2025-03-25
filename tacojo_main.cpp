#include "trans_ancestry_cojo.h"


void TransAncestryCOJO::read_files(string cojoFile1, string cojoFile2, string PLINK1, string PLINK2) 
{   
    list<string> commonSNP_temp;

    // read cojofiles and get rough common SNPs    
    c1.read_sumstat(cojoFile1);
    c2.read_sumstat(cojoFile2);
    
    // SNP already sorted
    set_intersection(c1.SNP.begin(), c1.SNP.end(), 
        c2.SNP.begin(), c2.SNP.end(), back_inserter(commonSNP_temp));
    
    // exclude rare SNP
    int i = 0;
    for (auto iter = commonSNP_temp.begin(); iter != commonSNP_temp.end(); ) {
        if (abs(c1.freq[c1.SNP_index[*iter]]-0.5) > 0.49 && abs(c2.freq[c2.SNP_index[*iter]]-0.5) > 0.49)
            iter = commonSNP_temp.erase(iter);
        else {
            Cohort::sumstat_commonSNP_index.insert(Cohort::sumstat_commonSNP_index.end(), pair<string, int> (*iter, i));
            iter++;
            i++;
        }
    }

    // set bim file 1 as reference, compare and get common SNPs across 4 files
    Cohort::A1_ref.resize(i);
    Cohort::A2_ref.resize(i);
    Cohort::SNP_pos_ref.resize(i);

    c1.read_PLINK(PLINK1, true);
    c2.read_PLINK(PLINK2, false);
    
    if (commonSNP_temp.size() != c1.included_SNP_PLINK.size() || 
        commonSNP_temp.size() != c2.included_SNP_PLINK.size()) {

        list<string>().swap(commonSNP_temp);
        
        // SNP already sorted
        set_intersection(c1.included_SNP_PLINK.begin(), c1.included_SNP_PLINK.end(), 
            c2.included_SNP_PLINK.begin(), c2.included_SNP_PLINK.end(), back_inserter(commonSNP_temp));
        
        vector<int> SNP_pos_ref_temp;

        for (auto iter = commonSNP_temp.begin(); iter != commonSNP_temp.end(); iter++)
            SNP_pos_ref_temp.push_back(Cohort::SNP_pos_ref[Cohort::sumstat_commonSNP_index[*iter]]);
        
        SNP_pos_ref_temp.swap(Cohort::SNP_pos_ref);
        vector<int>().swap(SNP_pos_ref_temp);        
    }

    if (commonSNP_temp.size() == 0)
        LOGGER.e(0, "Input data has no common SNPs.");

    LOGGER.i(0, "common SNPs included for two cohorts", to_string(commonSNP_temp.size()));
    LOGGER.i(0, "individuals in Cohort 1", to_string(c1.indi_num));
    LOGGER.i(0, "individuals in Cohort 2", to_string(c2.indi_num));

    move(commonSNP_temp.begin(), commonSNP_temp.end(), back_inserter(Cohort::commonSNP));
    list<string>().swap(commonSNP_temp);

    // initialize sumstat and X matrices, calculate Vp
    c1.generate_sumstat_and_X();
    c2.generate_sumstat_and_X();
    cout << "Vp: " << c1.Vp << " " << c2.Vp << endl << endl;
}


void TransAncestryCOJO::initialize_hyperparameters()
{   
    threshold = 5e-8;
    window_size = INT_MAX; // 1500000;

    colinear_threshold = 0.8;
    colinear_threshold_sqrt = sqrt(colinear_threshold);
    iter_colinear_threshold = 1 / (1 - colinear_threshold);
    
    R2_incremental_threshold = 0.0;
    R2_incremental_threshold_backwards = -0.5;
    max_iter_num = 1000;  
}
    

void TransAncestryCOJO::append_row(ArrayXXd &matrix, const ArrayXXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void TransAncestryCOJO::append_row(MatrixXd &matrix, const MatrixXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void TransAncestryCOJO::append_column(MatrixXd &matrix, const MatrixXd &vector)
{
    int numCols = matrix.cols();
    matrix.conservativeResize(NoChange, numCols+1);
    matrix.col(numCols) = vector;
}


void TransAncestryCOJO::remove_row(ArrayXXd &matrix, int index)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void TransAncestryCOJO::remove_row(MatrixXd &matrix, int index)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void TransAncestryCOJO::remove_column(MatrixXd &matrix, int index)
{   
    // -1 indicates the last column
    int numRows = matrix.rows(), numCols = matrix.cols()-1;

    if (index != -1 && index < numCols)
        matrix.middleCols(index, numCols-index) = matrix.rightCols(numCols-index).eval();

    matrix.conservativeResize(NoChange, numCols);
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


void TransAncestryCOJO::calc_inner_product(const vector<int> &index_list, int single_index) 
{   
    c1.r_temp_vec.resize(index_list.size());
    c2.r_temp_vec.resize(index_list.size());

    # pragma omp parallel for 
    for (int i = 0; i < index_list.size(); i++) {
        if (abs(Cohort::SNP_pos_ref[index_list[i]] - Cohort::SNP_pos_ref[single_index]) <= window_size) {
            c1.r_temp_vec(i) = (c1.X.col(single_index).transpose() * c1.X.col(index_list[i])).value() / (c1.indi_num - 1);
            c2.r_temp_vec(i) = (c2.X.col(single_index).transpose() * c2.X.col(index_list[i])).value() / (c2.indi_num - 1);
        }
        else {
            c1.r_temp_vec(i) = 0;
            c2.r_temp_vec(i) = 0;
        }
    }
}


void TransAncestryCOJO::initialize_matrices() 
{   
    candidate_SNP.push_back(max_SNP_index);
    screened_SNP.resize(Cohort::commonSNP.size());
    iota(screened_SNP.begin(), screened_SNP.end(), 0);

    calc_inner_product(screened_SNP, max_SNP_index);

    c1.sumstat_candidate = c1.sumstat.row(max_SNP_index);
    c2.sumstat_candidate = c2.sumstat.row(max_SNP_index);

    c1.sumstat_screened = c1.sumstat;
    c2.sumstat_screened = c2.sumstat;

    c1.r = c1.r_temp_vec;
    c2.r = c2.r_temp_vec;

    c1.R_inv_pre = MatrixXd::Identity(1,1);
    c2.R_inv_pre = MatrixXd::Identity(1,1);
}


void TransAncestryCOJO::remove_new_colinear_SNP() 
{   
    for (int i = screened_SNP.size()-1, last_col = candidate_SNP.size()-1; i >= 0; i--) {
        if (abs(c1.r(i, last_col)) >= colinear_threshold_sqrt || abs(c2.r(i, last_col)) >= colinear_threshold_sqrt) {
            remove_row(c1.sumstat_screened, i);
            remove_row(c2.sumstat_screened, i);
            remove_row(c1.r, i);
            remove_row(c2.r, i);
            
            if (i != max_SNP_index)
                excluded_SNP.push_back(screened_SNP[i]);
            screened_SNP.erase(screened_SNP.begin()+i);
        }
    }
}   


void TransAncestryCOJO::refuse_max_SNP_as_candidate() 
{
    remove_row(c1.sumstat_candidate);
    remove_row(c2.sumstat_candidate);
    sumstat_merge(max_SNP_index, 2) = -1;
}


void TransAncestryCOJO::initialize_backward_selection() 
{   
    c1.sumstat_backward = c1.sumstat_candidate;
    c2.sumstat_backward = c2.sumstat_candidate;

    for (int i = candidate_SNP.size()-1; i >= 0; i--) {
        if (sumstat_new_model_joint(i, 3) > threshold) {
            remove_row(c1.sumstat_backward, i);
            remove_row(c2.sumstat_backward, i);
        }
    }

    vector<VectorXd> X_backward_cohort1, X_backward_cohort2;
    vector<int> SNP_pos_backward;
    int total_num = 0;
    
    for (int i = 0; i < candidate_SNP.size(); i++) {
        if (sumstat_new_model_joint(i, 3) <= threshold) {
            X_backward_cohort1.push_back(c1.X.col(candidate_SNP[i]));
            X_backward_cohort2.push_back(c2.X.col(candidate_SNP[i]));
            SNP_pos_backward.push_back(Cohort::SNP_pos_ref[candidate_SNP[i]]);
            total_num++;
        }
    }

    MatrixXd R_post_cohort1 = MatrixXd::Identity(total_num, total_num);
    MatrixXd R_post_cohort2 = MatrixXd::Identity(total_num, total_num);

    if (total_num > 1) {
        // # pragma omp parallel for 
        for (int i = 1; i < total_num; i++)
            for (int j = 0; j < i; j++)
                if (abs(SNP_pos_backward[i] - SNP_pos_backward[j]) <= window_size) {
                    R_post_cohort1(i, j) = (X_backward_cohort1[i].transpose() * X_backward_cohort1[j]).value() / (c1.indi_num-1);
                    R_post_cohort1(j, i) = R_post_cohort1(i, j);

                    R_post_cohort2(i, j) = (X_backward_cohort2[i].transpose() * X_backward_cohort2[j]).value() / (c2.indi_num-1);
                    R_post_cohort2(j, i) = R_post_cohort2(i, j);
                } 
    }

    vector<VectorXd>().swap(X_backward_cohort1);
    vector<VectorXd>().swap(X_backward_cohort2);
        
    c1.R_inv_post.noalias() = R_post_cohort1.inverse();
    c2.R_inv_post.noalias() = R_post_cohort2.inverse();
}


void TransAncestryCOJO::adjust_SNP_according_to_backward_selection() 
{   
    remove_row(c1.sumstat_screened, max_SNP_index);
    remove_row(c2.sumstat_screened, max_SNP_index);

    remove_row(c1.r, max_SNP_index);
    remove_row(c2.r, max_SNP_index);

    screened_SNP.erase(screened_SNP.begin()+max_SNP_index);

    for (int i = candidate_SNP.size()-1; i >= 0; i--) {
        if (sumstat_new_model_joint(i, 3) > threshold) {
            remove_row(c1.sumstat_candidate, i);
            remove_row(c2.sumstat_candidate, i);
            remove_column(c1.r, i);
            remove_column(c2.r, i);  

            LOGGER.w(1, "Previous candidate SNP removed", Cohort::commonSNP[candidate_SNP[i]]);
            backward_removed_SNP.push_back(candidate_SNP[i]);
            candidate_SNP.erase(candidate_SNP.begin()+i);
        }
    }
    
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {
        calc_inner_product(candidate_SNP, excluded_SNP[i]); 
        
        if (c1.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt && 
            c2.r_temp_vec.cwiseAbs().maxCoeff() < colinear_threshold_sqrt) {
            append_row(c1.sumstat_screened, c1.sumstat.row(excluded_SNP[i]));
            append_row(c2.sumstat_screened, c2.sumstat.row(excluded_SNP[i])); 
            append_row(c1.r, c1.r_temp_vec.transpose());
            append_row(c2.r, c2.r_temp_vec.transpose());
            
            screened_SNP.push_back(excluded_SNP[i]); 
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }
}


void TransAncestryCOJO::main_loop() 
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

    double previous_R2_cohort1 = 0.0, previous_R2_cohort2 = 0.0;
    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;

    initialize_matrices();
    remove_new_colinear_SNP();

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        cout << candidate_SNP.size() << " " << screened_SNP.size() << " " << backward_removed_SNP.size() << endl;

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
            
            string temp_SNP_name = Cohort::commonSNP[screened_SNP[max_SNP_index]];

            // calculate joint effects
            append_row(c1.sumstat_candidate, c1.sumstat_screened.row(max_SNP_index));
            append_row(c2.sumstat_candidate, c2.sumstat_screened.row(max_SNP_index));

            calc_inner_product(candidate_SNP, screened_SNP[max_SNP_index]);
            c1.calc_R_inv_fast();
            c2.calc_R_inv_fast();

            NA_flag = c1.calc_joint_effects(c1.sumstat_candidate, true, iter_colinear_threshold);       
            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }
            
            NA_flag = c2.calc_joint_effects(c2.sumstat_candidate, true, iter_colinear_threshold);
            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }

            inverse_var_meta(c1.beta, c2.beta, c1.beta_var, c2.beta_var, sumstat_new_model_joint);

            if ((c1.R2 < (1+R2_incremental_threshold) * previous_R2_cohort1) || 
                (c2.R2 < (1+R2_incremental_threshold) * previous_R2_cohort2)) {
                // LOGGER.w(1, "R2 increment unsatisfactory", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }

            // include new candidate SNP
            calc_inner_product(screened_SNP, screened_SNP[max_SNP_index]);

            candidate_SNP.push_back(screened_SNP[max_SNP_index]);
            append_column(c1.r, c1.r_temp_vec);
            append_column(c2.r, c2.r_temp_vec);
            
            if (sumstat_new_model_joint.col(3).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed", temp_SNP_name);

                int M = candidate_SNP.size();
                // cout << "Added R vector cohort 1: " << c1.r.row(max_SNP_index) << endl;
                // cout << "Added R vector cohort 2: " << c2.r.row(max_SNP_index) << endl;
                remove_new_colinear_SNP();
                cout << "Added diagonal value cohort 1: " << c1.R_inv_post(M-1, M-1) << endl;
                cout << "Added diagonal value cohort 2: " << c2.R_inv_post(M-1, M-1) << endl;
                cout << "Joint b: " << sumstat_new_model_joint(M-1, 0) << endl;
                cout << "Joint se: " << sqrt(sumstat_new_model_joint(M-1, 1)) << endl;
                cout << "Joint p-value: " << scientific << sumstat_new_model_joint(M-1, 3) << endl;
                cout << "Adjusted R2 for cohort 1: " << fixed << c1.R2 << endl;
                cout << "Adjusted R2 for cohort 2: " << c2.R2 << endl;
                break; 
            }

            // backward selection
            LOGGER << "Backward selection" << endl;

            initialize_backward_selection();
            c1.calc_joint_effects(c1.sumstat_backward, false);
            c2.calc_joint_effects(c2.sumstat_backward, false);

            cout << "Adjusted R2 after SNP elimination for cohort 1: " << fixed << c1.R2 << endl;
            cout << "Adjusted R2 after SNP elimination for cohort 2: " << c2.R2 << endl;

            if ((c1.R2 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort1) || 
                (c2.R2 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort2)) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", temp_SNP_name);
                candidate_SNP.pop_back();
                remove_column(c1.r);
                remove_column(c2.r);
                refuse_max_SNP_as_candidate();
                continue;
            }  
            
            LOGGER.d(0, "Backward selection succeeded", temp_SNP_name);
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

        cout << "iter " << ++iter_num << " finished" << endl;
        cout << "--------------------------------" << endl;
    }

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);
}


void TransAncestryCOJO::save_results(string filepath) 
{
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    vector<string> headers = {"SNP", "A1", "A2", 
        "freq.x", "b.x", "se.x", "p.x", "N.x", "D.x", "V.x", "bJ.x", "seJ.x", "zJ.x", "pJ.x", 
        "freq.y", "b.y", "se.y", "p.y", "N.y", "D.y", "V.y", "bJ.y", "seJ.y", "zJ.y", "pJ.y", 
        "seJ.ma", "bJ.ma", "zJ", "pJ"};
        
    for (auto iter = headers.begin(); iter != headers.end(); iter++)
        jmaCOJO << *iter << "\t";
        
    jmaCOJO << endl;

    // ItemBim temp_item;
    string SNP_name;
    int index = 0;
    map<string, int> SNP_order_pair;
 
    // Inserting element in pair vector
    // to keep track of previous indexes
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

    for (auto iter = SNP_order_pair.begin(); iter != SNP_order_pair.end(); iter++) {
        SNP_name = iter->first;
        index = iter->second;
        jmaCOJO << SNP_name << "\t" << Cohort::A1_ref[Cohort::sumstat_commonSNP_index[SNP_name]] << "\t" 
                << Cohort::A2_ref[Cohort::sumstat_commonSNP_index[SNP_name]] << "\t";

        // cols 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
        temp_row = c1.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t" << temp_row(6) << "\t" << temp_row(5) << "\t";
        jmaCOJO << bJx(index) << "\t" << seJx(index) << "\t" << zJx(index) << "\t" << pJx(index) << "\t";
            
        temp_row = c2.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t" << temp_row(6) << "\t" << temp_row(5) << "\t";
        jmaCOJO << bJy(index) << "\t" << seJy(index) << "\t" << zJy(index) << "\t" << pJy(index) << "\t";
        
        jmaCOJO << bJma(index) << "\t" << seJma(index) << "\t" << zJma(index) << "\t" << pJma(index) << endl;
    }

    LOGGER.i(0, "Results saved into [" + filepath + "]", "Finished");
    jmaCOJO.close();
}