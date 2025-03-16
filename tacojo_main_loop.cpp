#include "trans_ancestry_cojo.h"


void taCOJO::initialize_hyperparameters()
{   
    threshold = 5e-8;

    colinearity_threshold = 0.8;
    colinearity_threshold_sqrt = sqrt(colinearity_threshold);
    iter_colinearity_threshold = 1 / (1 - colinearity_threshold);
    
    R2_incremental_threshold = 0.0;
    R2_incremental_threshold_backwards = -0.5;
    max_iter_num = 10000;  
}
    

void taCOJO::inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
    const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge) 
{
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    merge.resize(b_cohort1.size(), 4);

    merge.col(0) = (b_cohort1 * se2_cohort2 + b_cohort2 * se2_cohort1) / (se2_cohort1 + se2_cohort2);
    merge.col(1) = se2_cohort1 * se2_cohort2 / (se2_cohort1 + se2_cohort2);
    merge.col(2) = abs(merge.col(0) / sqrt(merge.col(1)));
    merge.col(3) = erfc(merge.col(2)/sqrt(2));
}


void taCOJO::append_row(ArrayXXd &matrix, const ArrayXXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void taCOJO::append_row(MatrixXd &matrix, const MatrixXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void taCOJO::append_column(MatrixXd &matrix, const MatrixXd &vector)
{
    int numCols = matrix.cols();
    matrix.conservativeResize(NoChange, numCols+1);
    matrix.col(numCols) = vector;
}


void taCOJO::remove_row(ArrayXXd &matrix, int index=-1)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void taCOJO::remove_row(MatrixXd &matrix, int index=-1)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void taCOJO::remove_column(MatrixXd &matrix, int index=-1)
{   
    // -1 indicates the last column
    int numRows = matrix.rows(), numCols = matrix.cols()-1;

    if (index != -1 && index < numCols)
        matrix.middleCols(index, numCols-index) = matrix.rightCols(numCols-index).eval();

    matrix.conservativeResize(NoChange, numCols);
}


void taCOJO::initialize_matrices() 
{   
    candidate_SNP.push_back(max_SNP_index);
    screened_SNP.resize(commonSNP_num);
    iota(screened_SNP.begin(), screened_SNP.end(), 0);

    sumstat1_candidate = sumstat1.row(max_SNP_index);
    sumstat2_candidate = sumstat2.row(max_SNP_index);
    X1_candidate = X1.col(max_SNP_index);
    X2_candidate = X2.col(max_SNP_index);

    sumstat1_screened = sumstat1;
    sumstat2_screened = sumstat2;
    X1_screened = X1;
    X2_screened = X2;

    r1 = X1_candidate.transpose() * X1_screened / (indi_num1 - 1);
    r2 = X2_candidate.transpose() * X2_screened / (indi_num2 - 1);
}


void taCOJO::remove_new_colinear_SNP() 
{   
    for (int i = screened_SNP.size()-1, last_row = candidate_SNP.size()-1; i >= 0; i--) {
        if (abs(r1(last_row, i)) >= colinearity_threshold_sqrt || abs(r2(last_row, i)) >= colinearity_threshold_sqrt) {
            remove_row(sumstat1_screened, i);
            remove_row(sumstat2_screened, i);
            remove_column(r1, i);
            remove_column(r2, i);
            remove_column(X1_screened, i);
            remove_column(X2_screened, i);
            
            if (i != max_SNP_index)
                excluded_SNP.push_back(screened_SNP[i]);
            screened_SNP.erase(screened_SNP.begin()+i);
        }
    }
}   


void taCOJO::calc_conditional_effects(ArrayXd &conditional_beta1, ArrayXd &conditional_beta2) 
{   
    MatrixXd temp1, temp2;
    ArrayXd temp3;

    temp1 = sumstat1_candidate.col(0) * sqrt(sumstat1_candidate.col(5));
    temp2.noalias() = R1_inv_pre * temp1;
    temp3 = r1.transpose() * temp2;
    conditional_beta1 = sumstat1_screened.col(0) - temp3 / sqrt(sumstat1_screened.col(5));
    
    temp1 = sumstat2_candidate.col(0) * sqrt(sumstat2_candidate.col(5));
    temp2.noalias() = R2_inv_pre * temp1;
    temp3 = r2.transpose() * temp2;
    conditional_beta2 = sumstat2_screened.col(0) - temp3 / sqrt(sumstat2_screened.col(5));
}


void taCOJO::fast_inv(const MatrixXd &R_inv_pre, const VectorXd &new_column, MatrixXd &R_inv_post) {
    double temp_element = 1 / (1 - new_column.transpose() * R_inv_pre * new_column);
    VectorXd temp_vector = R_inv_pre * new_column;

    int dim = R_inv_pre.rows();
    R_inv_post.resize(dim+1, dim+1);

    R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
    R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
    R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
    R_inv_post(dim, dim) = temp_element;
}


bool taCOJO::calc_joint_effects(const ArrayXXd &sumstat_candidate, const MatrixXd &R_inv_post, 
    double Vp, ArrayXd &beta, ArrayXd &beta_var, double &R2, bool flag) 
{       
    if (flag && ((abs(R_inv_post.minCoeff()) > iter_colinearity_threshold) || (abs(R_inv_post.maxCoeff()) > iter_colinearity_threshold)))
        return true;

    VectorXd temp1 = sqrt(sumstat_candidate.col(5)) * sumstat_candidate.col(0);
    ArrayXd temp2 = R_inv_post * temp1;
    beta = temp2 / sqrt(sumstat_candidate.col(5));
    
    double sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
                
    beta_var = sigma_J_squared * R_inv_post.diagonal().array() / sumstat_candidate.col(6) ;

    if (flag && beta_var.minCoeff() <= 0)
        return true;

    double Neff = median(sumstat_candidate.col(4));
    int M = sumstat_candidate.rows();
    R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    return false;
}


void taCOJO::refuse_max_SNP_as_candidate() 
{
    // cohort 1
    remove_row(sumstat1_candidate);

    // cohort 2
    if (sumstat1_candidate.rows() != sumstat2_candidate.rows())
        remove_row(sumstat2_candidate);

    sumstat_merge(max_SNP_index, 2) = -1;
}


void taCOJO::initialize_backward_selection() 
{   
    sumstat1_new_model = sumstat1_candidate;
    sumstat2_new_model = sumstat2_candidate;
    X1_new_model = X1_candidate;
    X2_new_model = X2_candidate;

    for (int i = candidate_SNP.size()-1; i >= 0; i--) {
        if (new_model_joint(i, 3) > threshold) {
            remove_row(sumstat1_new_model, i);
            remove_row(sumstat2_new_model, i);
            remove_column(X1_new_model, i);
            remove_column(X2_new_model, i);
        }
    }

    R1_inv_post.noalias() = (X1_new_model.transpose() * X1_new_model).inverse() * (indi_num1 - 1);
    R2_inv_post.noalias() = (X2_new_model.transpose() * X2_new_model).inverse() * (indi_num2 - 1);
}


void taCOJO::adjust_SNP_according_to_backward_selection() 
{   
    remove_row(sumstat1_screened, max_SNP_index);
    remove_row(sumstat2_screened, max_SNP_index);
    remove_column(r1, max_SNP_index);
    remove_column(r2, max_SNP_index);
    remove_column(X1_screened, max_SNP_index);
    remove_column(X2_screened, max_SNP_index);
    screened_SNP.erase(screened_SNP.begin()+max_SNP_index);

    for (int i = candidate_SNP.size()-1; i >= 0; i--) {
        if (new_model_joint(i, 3) > threshold) {
            remove_row(sumstat1_candidate, i);
            remove_row(sumstat2_candidate, i);
            remove_row(r1, i);
            remove_row(r2, i);
            remove_column(X1_candidate, i);
            remove_column(X2_candidate, i);    

            LOGGER.w(1, "Previous candidate SNP removed", commonSNP[candidate_SNP[i]]);
            backward_removed_SNP.push_back(candidate_SNP[i]);
            candidate_SNP.erase(candidate_SNP.begin()+i);
        }
    }
    
    VectorXd r1_temp_vec, r2_temp_vec;
    int temp_index;
    
    for (int i = excluded_SNP.size()-1; i >= 0; i--) {
        temp_index = excluded_SNP[i];
        r1_temp_vec = X1.col(temp_index).transpose() * X1_candidate / (indi_num1 - 1);
        r2_temp_vec = X2.col(temp_index).transpose() * X2_candidate / (indi_num2 - 1);
        
        if (r1_temp_vec.cwiseAbs().maxCoeff() < colinearity_threshold_sqrt && 
            r2_temp_vec.cwiseAbs().maxCoeff() < colinearity_threshold_sqrt) {
            append_row(sumstat1_screened, sumstat1.row(temp_index));
            append_row(sumstat2_screened, sumstat2.row(temp_index));  
            append_column(X1_screened, X1.col(temp_index));
            append_column(X2_screened, X2.col(temp_index)); 
            append_column(r1, r1_temp_vec);
            append_column(r2, r2_temp_vec);

            screened_SNP.push_back(temp_index);
            excluded_SNP.erase(excluded_SNP.begin()+i);
        }
    }
}


void taCOJO::save_results(string filepath) 
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
        SNP_order_pair.insert(make_pair(commonSNP[*iter], index));
    
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
        jmaCOJO << SNP_name << "\t" << A1_ref[commonSNP_index[SNP_name]] << "\t" << A2_ref[commonSNP_index[SNP_name]] << "\t";

        // cols 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
        temp_row = sumstat1_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t" << temp_row(6) << "\t" << temp_row(5) << "\t";
        jmaCOJO << bJx(index) << "\t" << seJx(index) << "\t" << zJx(index) << "\t" << pJx(index) << "\t";
            
        temp_row = sumstat2_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t" << temp_row(6) << "\t" << temp_row(5) << "\t";
        jmaCOJO << bJy(index) << "\t" << seJy(index) << "\t" << zJy(index) << "\t" << pJy(index) << "\t";
        
        jmaCOJO << bJma(index) << "\t" << seJma(index) << "\t" << zJma(index) << "\t" << pJma(index) << endl;
    }

    LOGGER.i(0, "Results saved into [" + filepath + "]", "Finished");
    jmaCOJO.close();
}


void taCOJO::main_loop() 
{   
    output_b_cohort1 = sumstat1.col(0);
    output_b_cohort2 = sumstat2.col(0);
    output_se2_cohort1 = sumstat1.col(1);
    output_se2_cohort2 = sumstat2.col(1);
    R1_inv_pre = MatrixXd::Identity(1,1);
    R2_inv_pre = MatrixXd::Identity(1,1);

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);

    sumstat_merge.col(2).maxCoeff(&max_SNP_index);
    if (sumstat_merge(max_SNP_index, 3) > threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    LOGGER.i(0, "First SNP", commonSNP[max_SNP_index]);
    cout << "--------------------------------" << endl;

    ArrayXd conditional_beta1, conditional_beta2;
    ArrayXd beta1, beta2, beta_var1, beta_var2;
    double R2_cohort1, R2_cohort2;
    double previous_R2_cohort1 = 0.0, previous_R2_cohort2 = 0.0;

    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;

    initialize_matrices();
    remove_new_colinear_SNP();

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        cout << candidate_SNP.size() << " " << screened_SNP.size() << " " << backward_removed_SNP.size() << endl;
        // calculate conditional effects
        calc_conditional_effects(conditional_beta1, conditional_beta2);
        inverse_var_meta(conditional_beta1, conditional_beta2, sumstat1_screened.col(1), sumstat2_screened.col(1), sumstat_merge);

        while (true) {
            // select maximal SNP
            double max_Zabs = sumstat_merge.col(2).maxCoeff(&max_SNP_index);
            if (max_Zabs < 0 || sumstat_merge(max_SNP_index, 3) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }

            string temp_SNP_name = commonSNP[screened_SNP[max_SNP_index]];

            // calculate joint effects
            append_row(sumstat1_candidate, sumstat1_screened.row(max_SNP_index));
            append_row(sumstat2_candidate, sumstat2_screened.row(max_SNP_index));
        
            fast_inv(R1_inv_pre, X1_screened.col(max_SNP_index).transpose() * X1_candidate / (indi_num1 - 1), R1_inv_post);
            fast_inv(R2_inv_pre, X2_screened.col(max_SNP_index).transpose() * X2_candidate / (indi_num2 - 1), R2_inv_post);
            
            NA_flag = calc_joint_effects(sumstat1_candidate, R1_inv_post, Vp1, beta1, beta_var1, R2_cohort1, true);            
            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }
            
            NA_flag = calc_joint_effects(sumstat2_candidate, R2_inv_post, Vp2, beta2, beta_var2, R2_cohort2, true);
            if (NA_flag) {
                // LOGGER.w(1, "NA produced, potentially due to colinearity", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }

            inverse_var_meta(beta1, beta2, beta_var1, beta_var2, new_model_joint);

            if ((R2_cohort1 < (1+R2_incremental_threshold) * previous_R2_cohort1) || 
                (R2_cohort2 < (1+R2_incremental_threshold) * previous_R2_cohort2)) {
                // LOGGER.w(1, "R2 increment unsatisfactory", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }

            // include new candidate SNP
            candidate_SNP.push_back(screened_SNP[max_SNP_index]);
            append_column(X1_candidate, X1_screened.col(max_SNP_index));
            append_column(X2_candidate, X2_screened.col(max_SNP_index));
            append_row(r1, X1_screened.col(max_SNP_index).transpose() * X1_screened / (indi_num1 - 1));
            append_row(r2, X2_screened.col(max_SNP_index).transpose() * X2_screened / (indi_num2 - 1));

            if (new_model_joint.col(3).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed", temp_SNP_name);
                remove_new_colinear_SNP();

                // cout << "Added R vector cohort 1: " << r1.row(max_SNP_index) << endl;
                // cout << "Added R vector cohort 2: " << r2.row(max_SNP_index) << endl;
                int M = candidate_SNP.size();
                cout << "Added diagnoal value cohort 1: " << R1_inv_post(M-1, M-1) << endl;
                cout << "Added diagnoal value cohort 2: " << R2_inv_post(M-1, M-1) << endl;
                cout << "Joint b: " << new_model_joint(M-1, 0) << endl;
                cout << "Joint se: " << sqrt(new_model_joint(M-1, 1)) << endl;
                cout << "Joint p-value: " << scientific << new_model_joint(M-1, 3) << endl;
                cout << "Adjusted R2 for cohort 1: " << fixed << R2_cohort1 << endl;
                cout << "Adjusted R2 for cohort 2: " << R2_cohort2 << endl;
                break; 
            }

            // backward selection
            initialize_backward_selection();
            calc_joint_effects(sumstat1_new_model, R1_inv_post, Vp1, beta1, beta_var1, R2_cohort1, false);
            calc_joint_effects(sumstat2_new_model, R2_inv_post, Vp2, beta2, beta_var2, R2_cohort2, false);

            cout << "Adjusted R2 after SNP elimination for cohort 1: " << fixed << R2_cohort1 << endl;
            cout << "Adjusted R2 after SNP elimination for cohort 2: " << R2_cohort2 << endl;

            if ((R2_cohort1 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort1) || 
                (R2_cohort2 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort2)) {
                LOGGER.w(1, "Backward selection, adjusted R2 lower than threshold", temp_SNP_name);
                candidate_SNP.pop_back();
                remove_column(X1_candidate);
                remove_column(X2_candidate);
                remove_row(r1);
                remove_row(r2);           
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
            previous_R2_cohort1 = R2_cohort1;
            previous_R2_cohort2 = R2_cohort2;
            R1_inv_pre = R1_inv_post;
            R2_inv_pre = R2_inv_post;

            output_b_cohort1 = beta1;
            output_b_cohort2 = beta2;
            output_se2_cohort1 = beta_var1;
            output_se2_cohort2 = beta_var2; 
        }

        cout << "iter " << ++iter_num << " finished" << endl;
        cout << "--------------------------------" << endl;
    }

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);
}
