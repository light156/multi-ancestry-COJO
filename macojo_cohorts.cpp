#include "macojo.h"


void Cohort::get_vector_from_bed_matrix(int index, ArrayXd &vec)
{
    vec.setZero(indi_num);
    uint64_t X_start_index = uint64_t(index) * indi_num;
    uint64_t start_word = X_start_index >> 6;
    uint64_t start_bit  = X_start_index & 63;

    double SNP_avg = X_avg[index], SNP_std = X_std[index];
    int i = 0;

    while (i < indi_num) {
        uint64_t wordA1 = X_A1.data_[start_word];
        uint64_t wordA2 = X_A2.data_[start_word];

        for (uint64_t b = start_bit; b < 64 && i < indi_num; b++, i++) {
            bool A1 = (wordA1 >> b) & 1ULL;
            bool A2 = (wordA2 >> b) & 1ULL;

            if (if_gcta_COJO || if_fill_NA)
                vec(i) = (!A1 || A2) ? (double(A1) + double(A2) - SNP_avg) / SNP_std : 0;
            else
                vec(i) = (!A1 || A2) ? double(A1) + double(A2) : -9;
        }

        start_word++;
        start_bit = 0;
    }
}


void Cohort::calc_inner_product_with_SNP_list(const vector<int> &SNP_list, const ArrayXXd &sumstat_temp, int single_index)
{      
    if (SNP_list.size() != sumstat_temp.rows()) {
        // this should not happen
        LOGGER.e(0, "Error: SNP list size not equal to sumstat_temp rows");
        return;
    }

    r_temp_vec.setZero(SNP_list.size());
    int temp_index = MACOJO::final_commonSNP_index[single_index];
    int temp_SNP_pos = MACOJO::SNP_pos_ref[temp_index];

    if (if_LD_mode) {
        for (int j = 0; j < SNP_list.size(); j++) {
            if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]] - temp_SNP_pos) <= window_size)
                r_temp_vec(j) = LD_matrix(MACOJO::final_commonSNP_index[SNP_list[j]], temp_index);
        }
    } else {
        get_vector_from_bed_matrix(temp_index, X_temp_vec);

        #pragma omp parallel for
        for (int j = 0; j < SNP_list.size(); j++) {
            if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]] - temp_SNP_pos) <= window_size) {
                ArrayXd vec_in_list;
                get_vector_from_bed_matrix(MACOJO::final_commonSNP_index[SNP_list[j]], vec_in_list);
                r_temp_vec(j) = calc_inner_product(X_temp_vec, vec_in_list, if_fill_NA);
            }       
        }
    }

    r_temp_vec_plain = r_temp_vec;

    if (if_gcta_COJO) {
        double N_candidate = sumstat(single_index, 4);
        double V_candidate = sumstat(single_index, 5);
        r_temp_vec = r_temp_vec.array() * sumstat_temp.col(4).min(N_candidate) * \
            sqrt(sumstat_temp.col(5)) * sqrt(V_candidate);
    }
}


void Cohort::calc_conditional_effects() 
{   
    MatrixXd temp1, temp2;
    ArrayXd temp3;

    if (if_gcta_COJO) {
        temp1 = sumstat_candidate.col(0) * sumstat_candidate.col(6);
        temp2.noalias() = R_inv_pre * temp1;
        temp3 = r * temp2;
        conditional_beta = sumstat_screened.col(0) - temp3 / sumstat_screened.col(6);
    } else {
        temp1 = sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5));
        temp2.noalias() = R_inv_pre * temp1;
        temp3 = r * temp2;
        conditional_beta = sumstat_screened.col(0) - temp3 / sqrt(sumstat_screened.col(5));
    }
}


void Cohort::calc_joint_effects(const ArrayXXd &s) 
{   
    if (if_gcta_COJO) {
        VectorXd temp1 = s.col(6) * s.col(0);
        beta = R_inv_post * temp1;
        beta_var = R_inv_post.diagonal().array() * Vp;
    } else {
        VectorXd temp1 = sqrt(s.col(5)) * s.col(0);
        ArrayXd temp2 = R_inv_post * temp1;
        beta = temp2 / sqrt(s.col(5));

        double sigma_J_squared = Vp - (s.col(0) * s.col(5) * beta).sum();
        beta_var = sigma_J_squared * R_inv_post.diagonal().array() / s.col(6);

        double Neff = median(s.col(4));
        int M = s.rows();
        R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    }
}


bool Cohort::calc_R_inv() {

    double lower_right = if_gcta_COJO ? sumstat_candidate(sumstat_candidate.rows()-1, 6) : 1;

    if (if_fast_inv) {
        VectorXd temp_vector = R_inv_pre * r_temp_vec;
        double temp_element = 1 / (lower_right - r_temp_vec.transpose() * temp_vector);
    
        int dim = R_inv_pre.rows();
        R_inv_post.setZero(dim+1, dim+1);
    
        R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
        R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
        R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
        R_inv_post(dim, dim) = temp_element;    
        return (abs(R_inv_post.minCoeff()) < iter_colinear_threshold) && (abs(R_inv_post.maxCoeff()) < iter_colinear_threshold);
    } else {
        int dim = R_pre.rows();
        R_post.setZero(dim+1, dim+1);
    
        R_post.block(0, 0, dim, dim) = R_pre;
        R_post.block(0, dim, dim, 1) = r_temp_vec;
        R_post.block(dim, 0, 1, dim) = r_temp_vec.transpose();
        R_post(dim, dim) = lower_right;
        /*
        LDLT<MatrixXd> ldlt_solver(R_post);
        VectorXd ldlt_D = ldlt_solver.vectorD();
        if (ldlt_D.minCoeff() <= 0 || sqrt(ldlt_D.maxCoeff() / ldlt_D.minCoeff()) > 30) 
            return false;
        */
        R_inv_post.noalias() = R_post.ldlt().solve(MatrixXd::Identity(dim+1, dim+1));
        return true;
    }
}


bool Cohort::calc_R_inv_from_SNP_list(const vector<int> &SNP_list, const ArrayXXd &sumstat_temp) 
{   
    if (SNP_list.size() != sumstat_temp.rows()) {
        // this should not happen
        LOGGER.e(0, "Error: SNP list size not equal to sumstat_temp rows");
        return false;
    }

    int total_num = SNP_list.size();
    R_post = MatrixXd::Identity(total_num, total_num);

    if (if_LD_mode) {
        for (int i = 0; i < total_num; i++) {
            for (int j = 0; j < i; j++) {
                if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[i]]] - 
                        MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]]) <= window_size)
                    R_post(i, j) = LD_matrix(MACOJO::final_commonSNP_index[SNP_list[i]], MACOJO::final_commonSNP_index[SNP_list[j]]);
            }
        }
    } else {
        vector<ArrayXd> X_new_model;
        
        for (int index : SNP_list) {
            get_vector_from_bed_matrix(MACOJO::final_commonSNP_index[index], X_temp_vec);
            X_new_model.push_back(X_temp_vec);
        }

        for (int i = 0; i < total_num; i++) {
            for (int j = 0; j < i; j++) {
                if (abs(MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[i]]] - 
                        MACOJO::SNP_pos_ref[MACOJO::final_commonSNP_index[SNP_list[j]]]) <= window_size)                
                    R_post(i, j) = calc_inner_product(X_new_model[i], X_new_model[j], if_fill_NA);
            }
        }
    } 
    
    for (int i = 0; i < total_num; i++) {
        if (if_gcta_COJO) 
            R_post(i, i) = sumstat_temp(i, 6);

        for (int j = 0; j < i; j++) {
            if (if_gcta_COJO) {
                R_post(i, j) = R_post(i, j) * min(sumstat_temp(i, 4), sumstat_temp(j, 4)) * \
                    sqrt(sumstat_temp(i, 5) * sumstat_temp(j, 5));
            }

            R_post(j, i) = R_post(i, j);
        }
    }

    R_inv_post.noalias() = R_post.ldlt().solve(MatrixXd::Identity(total_num, total_num));
    return (abs(R_inv_post.minCoeff()) < iter_colinear_threshold) && (abs(R_inv_post.maxCoeff()) < iter_colinear_threshold);
}
