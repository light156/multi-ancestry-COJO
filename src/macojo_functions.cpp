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


// make sure sumstat_candidate is prepared before calling this function
bool Cohort::calc_joint_effects(bool if_gcta_COJO)
{   
    if (if_gcta_COJO) {
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


bool Cohort::calc_R_inv_forward(int append_index)
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

    if (params.if_gcta_COJO) {
        R_inv_post_gcta = R_inv_pre_gcta - R_inv_pre_gcta.col(remove_index) * R_inv_pre_gcta.row(remove_index) / R_inv_pre_gcta(remove_index, remove_index);
        remove_row(R_inv_post_gcta, remove_index);
        remove_column(R_inv_post_gcta, remove_index);
    }
}


bool Cohort::calc_R_inv_from_SNP_list(const vector<int> &SNP_list, bool if_gcta_COJO, bool if_remove_NA) 
{   
    sumstat_candidate = sumstat(SNP_list, all);

    int total_num = SNP_list.size();
    R_post = MatrixXd::Identity(total_num, total_num);

    for (int i = 0; i < total_num; i++) {
        for (int j = 0; j < i; j++) {
            if (shared.chr_ref[SNP_list[i]] == shared.chr_ref[SNP_list[j]] &&
                abs(shared.SNP_pos_ref[SNP_list[i]] - shared.SNP_pos_ref[SNP_list[j]]) < params.window_size) {
                
                if (params.if_LD_mode)
                    R_post(i, j) = LD_matrix(SNP_list[i], SNP_list[j]);
                else
                    R_post(i, j) = genotype.calc_inner_product(SNP_list[i], SNP_list[j], if_remove_NA);

                R_post(j, i) = R_post(i, j);
            }
        }
    }
    
    R_inv_post = R_post.ldlt().solve(MatrixXd::Identity(total_num, total_num));

    if (if_gcta_COJO) {
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


void Cohort::update_r(const vector<int>& screened_SNP, int append_index) 
{  
    VectorXd r_temp_vec = VectorXd::Zero(shared.total_SNP_num);
    calc_inner_product_with_SNP_list(screened_SNP, append_index, r_temp_vec);
    append_column(r, r_temp_vec);

    if (params.if_gcta_COJO) {
        double N_new = sumstat(append_index, 4);
        double V_new = sumstat(append_index, 5);
        VectorXd r_temp_vec_gcta = r_temp_vec.array() * sumstat.col(4).min(N_new) * sqrt(sumstat.col(5)) * sqrt(V_new);
        append_column(r_gcta, r_temp_vec_gcta);
    }
}


void MACOJO::inverse_var_meta_conditional() 
{   
    bC = ArrayXd::Zero(shared.total_SNP_num);
    se2C = ArrayXd::Zero(shared.total_SNP_num);

    for (int n : current_list) {
        auto &c = cohorts[n]; 

        if (candidate_SNP.size() == 0) {
            bC += c.sumstat.col(0) / c.sumstat.col(1);
            se2C += 1 / c.sumstat.col(1);
        }
        else if (params.if_gcta_COJO) {
            VectorXd temp = c.R_inv_pre_gcta * (c.sumstat_candidate.col(0) * c.sumstat_candidate.col(6)).matrix();
            ArrayXd conditional_beta = c.sumstat.col(0) - (c.r_gcta * temp).array() / c.sumstat.col(6);

            bC += conditional_beta * c.sumstat.col(6) / c.Vp;
            se2C += c.sumstat.col(6) / c.Vp;
        }
        else {
            VectorXd temp = c.R_inv_pre * (c.sumstat_candidate.col(0) * sqrt(c.sumstat_candidate.col(5))).matrix();
            ArrayXd conditional_beta = c.sumstat.col(0) - (c.r * temp).array() / sqrt(c.sumstat.col(5));

            bC += conditional_beta / c.sumstat.col(1);
            se2C += 1 / c.sumstat.col(1);
        }
    }   
    
    se2C = 1 / se2C;
    bC = bC * se2C;
    abs_zC = abs(bC) / sqrt(se2C);

    abs_zC(removed_SNP) = -1;
    abs_zC(candidate_SNP) = -2;
    abs_zC(bad_SNP) = -3;
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


void MACOJO::output_cojo_cond(string savename) 
{   
    // this should never happen
    if (current_list.size() > 1)
        LOGGER.e(0, "Output conditional results is only available for single-cohort analysis");

    inverse_var_meta_conditional();
    auto& c = cohorts[current_list[0]];
    
    ofstream cmaCOJO(savename.c_str());
    if (!cmaCOJO) LOGGER.e(0, "cannot open the file [" + savename + "] to write.");

    cmaCOJO.precision(12);

    cmaCOJO << "Chr"
            << "\t" << "SNP" 
            << "\t" << "bp" 
            << "\t" << "A1" 
            << "\t" << "A2"
            << "\t" << "freq"
            << "\t" << "b"
            << "\t" << "se"
            << "\t" << "p"
            << "\t" << "n"
            << "\t" << "bC"
            << "\t" << "bC_se"
            << "\t" << "pC"
            << "\n";
    
    map<int, int> SNP_ref_order_pair;
    for (int index = 0; index < shared.total_SNP_num; index++)
        SNP_ref_order_pair.insert(make_pair(shared.SNP_pos_ref[index], index));

    for (auto iter = SNP_ref_order_pair.begin(); iter != SNP_ref_order_pair.end(); iter++) {
        int index = iter->second;
        
        // do not output bad SNPs that are not considered in the analysis
        // do not output candidate SNPs
        if (abs_zC(index) < -1.5) continue;

        cmaCOJO << shared.chr_ref[index]
                << "\t" << shared.SNP_ref[index]
                << "\t" << iter->first
                << "\t" << shared.A1_ref[index]
                << "\t" << shared.A2_ref[index];
        
        // cols 0:b, 1:se2, 2:p, 3:freq, 4:N 
        cmaCOJO << "\t" << c.sumstat(index, 3) 
                << "\t" << c.sumstat(index, 0) 
                << "\t" << sqrt(c.sumstat(index, 1)) 
                << "\t" << c.sumstat(index, 2)
                << "\t" << c.sumstat(index, 4);

        // output NA for removed SNPs
        if (abs_zC(index) < -0.5) 
            cmaCOJO << "\tNA\tNA\tNA\n";
        else 
            cmaCOJO << "\t" << bC(index)
                    << "\t" << sqrt(se2C(index))
                    << "\t" << erfc(abs_zC(index) / sqrt(2))
                    << "\n";
    }

    cmaCOJO.close();
    LOGGER.i(0, "Results saved into [" + savename + "]\n");
}



void MACOJO::output_cojo_joint(string savename) 
{
    ofstream jmaCOJO(savename.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + savename + "] to write.");

    jmaCOJO.precision(12);

    bool flag = true;
    
    for (int n : current_list) {
        if (!cohorts[n].calc_R_inv_from_SNP_list(candidate_SNP, params.effect_size_mode == "GCTA", params.effect_size_mode == "removeNA")) {
            LOGGER.w(0, "Colinearity exceeds threshold when calculating R inverse in Cohort " + to_string(n+1), "make sure this is what you want");
            flag = false;
        }

        cohorts[n].sumstat_candidate = cohorts[n].sumstat(candidate_SNP, all);

        if (!cohorts[n].calc_joint_effects(params.effect_size_mode == "GCTA")) {
            LOGGER.w(0, "Joint se too small in Cohort " + to_string(n+1), "make sure this is what you want");
            flag = false;
        }
    }

    if (!flag) savename += ".warning";

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
    }

    jmaCOJO << "\n";

    int index = 0;
    map<int, int> SNP_ref_order_pair;
    for (auto ref_index : candidate_SNP) {
        SNP_ref_order_pair.insert(make_pair(shared.SNP_pos_ref[ref_index], index));
        index++;
    }

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
                    << "\t" << c.beta(index) 
                    << "\t" << sqrt(c.beta_var(index)) 
                    << "\t" << erfc(abs(c.beta(index)) / sqrt(c.beta_var(index)) / sqrt(2));
        }

        if (current_list.size() > 1) {
            jmaCOJO << "\t" << bJ(index) 
                    << "\t" << sqrt(se2J(index))
                    << "\t" << erfc(abs_zJ(index) / sqrt(2));
        } 
        
        jmaCOJO << "\n";
    }

    jmaCOJO.close();
    LOGGER.i(0, "Results saved into [" + savename + "]\n");
}
