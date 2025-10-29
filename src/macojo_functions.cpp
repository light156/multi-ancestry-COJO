#include "macojo.h"


void Cohort::calc_inner_product_with_SNP_list(const vector<int> &SNP_list, int index, VectorXd &r_temp_vec)
{   
    r_temp_vec.setZero(shared.total_SNP_num);

    if (params.if_LD_mode) {
        for (size_t i = 0; i < SNP_list.size(); i++) {
            int sweep_index = SNP_list[i];
            if (shared.chr_ref[sweep_index] == shared.chr_ref[index] &&
                abs(shared.SNP_pos_ref[sweep_index] - shared.SNP_pos_ref[index]) < params.window_size)
                r_temp_vec(sweep_index) = LD_matrix(sweep_index, index);
        }
        return;
    }
    
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < SNP_list.size(); i++) {
        int sweep_index = SNP_list[i];
        if (shared.chr_ref[sweep_index] == shared.chr_ref[index] &&
            abs(shared.SNP_pos_ref[sweep_index] - shared.SNP_pos_ref[index]) < params.window_size) {
            r_temp_vec(sweep_index) = genotype.calc_inner_product(sweep_index, index, params.if_remove_NA);    
        }  
    }
}


void Cohort::calc_conditional_effects() 
{   
    if (params.if_gcta_COJO) {
        VectorXd temp = R_inv_pre_gcta * (sumstat_candidate.col(0) * sumstat_candidate.col(6)).matrix();
        conditional_beta = sumstat.col(0) - (r_gcta * temp).array() / sumstat.col(6);
    } else {
        VectorXd temp = R_inv_pre * (sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5))).matrix();
        conditional_beta = sumstat.col(0) - (r * temp).array() / sqrt(sumstat.col(5));
    }
}


bool Cohort::calc_joint_effects()
{   
    if (params.if_gcta_COJO) {
        beta = R_inv_post_gcta * (sumstat_candidate.col(0) * sumstat_candidate.col(6)).matrix();
        beta_var = R_inv_post_gcta.diagonal().array() * Vp;
    } else {
        ArrayXd temp = R_inv_post * (sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5))).matrix();
        beta = temp / sqrt(sumstat_candidate.col(5));

        double sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
        beta_var = sigma_J_squared * R_inv_post.diagonal().array() / sumstat_candidate.col(6);

        double Neff = median(sumstat_candidate.col(4));
        int M = sumstat_candidate.rows();
        R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    }

    return beta_var.minCoeff() > 1e-30;
}


void Cohort::save_temp_model() 
{   
    if (sumstat_candidate.rows() == 1) {
        output_b = sumstat_candidate.col(0);
        output_se2 = sumstat_candidate.col(1);
        R_inv_pre = MatrixXd::Identity(1,1);
        R_inv_pre_gcta = MatrixXd::Identity(1,1) / sumstat_candidate(0,6); // only used when if_gcta_COJO == true
        previous_R2 = -1; // only used for if_gcta_COJO == false;
    } else {
        output_b = beta;
        output_se2 = beta_var;
        R_inv_pre = R_inv_post;
        R_inv_pre_gcta = R_inv_post_gcta; // only used when if_gcta_COJO == true
        previous_R2 = R2; // only used for if_gcta_COJO == false
    }
}


bool Cohort::calc_R_inverse_forward(int append_index)
{   
    RowVectorXd r_temp_row = r.row(append_index);
    RowVectorXd temp_vector = r_temp_row * R_inv_pre;

    if (r_temp_row.cwiseAbs2().maxCoeff() > params.collinear && r_temp_row.dot(temp_vector) > params.collinear) 
        return false;

    double temp_element = 1.0 / (1.0 - r_temp_row.dot(temp_vector));

    int dim = R_inv_pre.rows();
    R_inv_post = R_inv_pre + temp_element * temp_vector.transpose() * temp_vector;
    R_inv_post.conservativeResize(dim + 1, dim + 1);
    R_inv_post.topRightCorner(dim, 1) = -temp_element * temp_vector.transpose();
    R_inv_post.bottomLeftCorner(1, dim) = -temp_element * temp_vector;
    R_inv_post(dim, dim) = temp_element;

    if (R_inv_post.cwiseAbs().maxCoeff() > params.iter_collinear_threshold)
        return false;

    if (params.if_gcta_COJO) {
        r_temp_row = r_gcta.row(append_index);
        temp_vector = r_temp_row * R_inv_pre_gcta;
        temp_element = 1.0 / (sumstat_candidate(sumstat_candidate.rows()-1, 6) - r_temp_row.dot(temp_vector));

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
void Cohort::calc_R_inverse_backward(int remove_index)
{
    int dim = R_inv_pre.rows();
    if (remove_index < 0 || remove_index >= dim)
        throw invalid_argument("calc_R_inverse_backward: remove_index out of range.");

    R_inv_post = R_inv_pre - R_inv_pre.col(remove_index) * R_inv_pre.row(remove_index) / R_inv_pre(remove_index, remove_index);
    remove_row(R_inv_post, remove_index);
    remove_column(R_inv_post, remove_index);

    if (params.if_gcta_COJO) {
        R_inv_post_gcta = R_inv_pre_gcta - R_inv_pre_gcta.col(remove_index) * R_inv_pre_gcta.row(remove_index) / R_inv_pre_gcta(remove_index, remove_index);
        remove_row(R_inv_post_gcta, remove_index);
        remove_column(R_inv_post_gcta, remove_index);
    }
}


bool Cohort::calc_R_inv_from_SNP_list(const vector<int> &SNP_list) 
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
                    R_post(i, j) = genotype.calc_inner_product(SNP_list[i], SNP_list[j], params.if_remove_NA);

                R_post(j, i) = R_post(i, j);
            }
        }
    }

    LDLT<MatrixXd> ldlt_solver(R_post);
    ArrayXd ldlt_D = ldlt_solver.vectorD().array();
    R_inv_post = ldlt_solver.solve(MatrixXd::Identity(total_num, total_num));

    if (params.if_gcta_COJO) {
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
    
    return ((ldlt_D.minCoeff() > 1e-8) && (sqrt(ldlt_D.maxCoeff() / ldlt_D.minCoeff()) < 30) 
            && (R_inv_post.diagonal().maxCoeff() < params.iter_collinear_threshold));
}


// must happen after sumstate_candidate is updated
void MACOJO::accept_SNP_as_candidate(int index) 
{   
    candidate_SNP.push_back(index);
    auto iter = find(screened_SNP.begin(), screened_SNP.end(), index);
    screened_SNP.erase(iter);

    for (int n : current_list) {
        auto &c = cohorts[n];
        
        VectorXd r_temp_vec = VectorXd::Zero(shared.total_SNP_num);
        c.calc_inner_product_with_SNP_list(screened_SNP, index, r_temp_vec);
        append_column(c.r, r_temp_vec);

        if (params.if_gcta_COJO) {
            double N_new = c.sumstat_candidate(c.sumstat_candidate.rows()-1, 4);
            double V_new = c.sumstat_candidate(c.sumstat_candidate.rows()-1, 5);
            VectorXd r_temp_vec_gcta = r_temp_vec.array() * c.sumstat.col(4).min(N_new) * \
                sqrt(c.sumstat.col(5)) * sqrt(V_new);
            append_column(c.r_gcta, r_temp_vec_gcta);
        }
    }
}


void MACOJO::remove_SNP_from_candidate(int candidate_index) 
{   
    removed_SNP.push_back(candidate_SNP[candidate_index]);
    candidate_SNP.erase(candidate_SNP.begin() + candidate_index);

    for (int n : current_list) {
        auto &c = cohorts[n];
        remove_row(c.sumstat_candidate, candidate_index);
        remove_column(c.r, candidate_index);

        if (params.if_gcta_COJO)
            remove_column(c.r_gcta, candidate_index);
    }
}


void MACOJO::inverse_var_meta_conditional() 
{   
    ArrayXd numer = ArrayXd::Zero(shared.total_SNP_num);
    ArrayXd denom = ArrayXd::Zero(shared.total_SNP_num);

    for (int n : current_list) {
        auto &c = cohorts[n]; 

        if (candidate_SNP.size() == 0) {
            numer += c.sumstat.col(0) / c.sumstat.col(1);
            denom += 1 / c.sumstat.col(1);
        }
        else if (params.if_gcta_COJO) {
            numer += c.conditional_beta * c.sumstat.col(6) / c.Vp;
            denom += c.sumstat.col(6) / c.Vp;
        }
        else {
            numer += c.conditional_beta / c.sumstat.col(1);
            denom += 1 / c.sumstat.col(1);
        }
    }   

    abs_zC = abs(numer) / sqrt(denom);
    abs_zC(bad_SNP) = -1;
    abs_zC(candidate_SNP) = -1;
    abs_zC(removed_SNP) = -1;
}


void MACOJO::inverse_var_meta_joint() 
{
    if (current_list.size() == 1) {
        auto &c = cohorts[current_list[0]];
        bJ = c.beta;
        se2J = c.beta_var;
    } else {
        bJ.setZero(cohorts[current_list[0]].beta.rows());
        se2J.setZero(cohorts[current_list[0]].beta.rows());

        for (int n : current_list) {
            auto &c = cohorts[n];
            bJ += c.beta / c.beta_var;
            se2J += 1 / c.beta_var;
        }

        se2J = 1 / se2J;
        bJ = bJ * se2J;
    }   

    abs_zJ = abs(bJ) / sqrt(se2J);
}


void MACOJO::output_results_to_file(string savename) 
{
    ofstream jmaCOJO(savename.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + savename + "] to write.");

    ArrayXd bJma, se2Jma;

    jmaCOJO << "Chr"
            << "\t" << "SNP" 
            << "\t" << "bp" 
            << "\t" << "A1" 
            << "\t" << "A2";

    if (current_list.size() == 1) {
        jmaCOJO << "\t" << "freq"
                << "\t" << "b"
                << "\t" << "se"
                << "\t" << "p"
                << "\t" << "n"
                << "\t" << "bJ"
                << "\t" << "bJ_se"
                << "\t" << "pJ";
    } else {
        for (int n : current_list) {
            jmaCOJO << "\t" << "freq." << to_string(n+1) 
                    << "\t" << "b." << to_string(n+1) 
                    << "\t" << "se." << to_string(n+1) 
                    << "\t" << "p." << to_string(n+1) 
                    << "\t" << "n." << to_string(n+1) 
                    << "\t" << "bJ." << to_string(n+1) 
                    << "\t" << "bJ_se." << to_string(n+1) 
                    << "\t" << "pJ." << to_string(n+1);
        }
        
        jmaCOJO << "\t" << "bJ.ma"  
                << "\t" << "bJ_se.ma" 
                << "\t" << "pJ.ma";

        bJma.setZero(candidate_SNP.size());
        se2Jma.setZero(candidate_SNP.size());

        for (int n : current_list) {
            auto &c = cohorts[n];
            bJma += c.output_b / c.output_se2;
            se2Jma += 1 / c.output_se2;
        }

        se2Jma = 1 / se2Jma; 
        bJma = bJma * se2Jma;
    }
    
    jmaCOJO << "\n";

    int index = 0;
    map<int, int> SNP_ref_order_pair;
    for (auto ref_index : candidate_SNP) {
        SNP_ref_order_pair.insert(make_pair(shared.SNP_pos_ref[ref_index], index));
        index++;
    }
    
    jmaCOJO.precision(12);

    for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) {
        index = iter->second;
        int ref_index = candidate_SNP[index];
        jmaCOJO << shared.chr_ref[ref_index]
                << "\t" << shared.SNP_ref[ref_index] 
                << "\t" << iter->first 
                << "\t" << shared.A1_ref[ref_index] 
                << "\t" << shared.A2_ref[ref_index];

        for (int n : current_list) {
            auto &c = cohorts[n];
            
            // cols 0:b, 1:se2, 2:p, 3:freq, 4:N, 5: bJ, 6:seJ, 7:pJ
            jmaCOJO << "\t" << c.sumstat_candidate(index, 3) 
                    << "\t" << c.sumstat_candidate(index, 0) 
                    << "\t" << sqrt(c.sumstat_candidate(index, 1)) 
                    << "\t" << c.sumstat_candidate(index, 2) 
                    << "\t" << c.sumstat_candidate(index, 4) 
                    << "\t" << c.output_b(index) 
                    << "\t" << sqrt(c.output_se2(index)) 
                    << "\t" << erfc(abs(c.output_b(index)) / sqrt(c.output_se2(index)) / sqrt(2));
        }

        if (current_list.size() > 1) {
            jmaCOJO << "\t" << bJma(index) 
                    << "\t" << sqrt(se2Jma(index))
                    << "\t" << erfc(abs(bJma(index)) / sqrt(se2Jma(index)) / sqrt(2));
        } 
        
        jmaCOJO << "\n";
    }

    jmaCOJO.clear();
    jmaCOJO.close();
    LOGGER.i(0, "Results saved into [" + savename + "]\n");
}
