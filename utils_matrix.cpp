#include "utils_matrix.h"


double calc_inner_product(const ArrayXd &vec1, const ArrayXd &vec2, bool if_keep_NA)
{
    if (!if_keep_NA)
        return (vec1 * vec2).sum();

    double s1=0, s2=0, s12=0, s11=0, s22=0, n=0;
    int indi_num = vec1.size();

    for (int k = 0; k < indi_num; k++) {
        if (vec1(k) > -5 && vec2(k) > -5) {
            s1 += vec1(k);
            s2 += vec2(k);
            s12 += vec1(k) * vec2(k);
            s11 += vec1(k) * vec1(k);
            s22 += vec2(k) * vec2(k);
            n++;
        }
    }

    return (n*s12 - s1*s2) / sqrt((n*s11 - s1*s1)*(n*s22 - s2*s2));
}


/*
bool calc_R_inverse_fast(
    const MatrixXd& R_inv_pre,
    const VectorXd& r_temp_vec,
    double lower_right_corner,
    double iter_colinear_threshold,
    MatrixXd& R_inv_post)
{   
    VectorXd temp_vector = R_inv_pre * r_temp_vec;
    double temp_element = 1.0 / (lower_right_corner - r_temp_vec.dot(temp_vector));

    int dim = R_inv_pre.rows();
    R_inv_post.setZero(dim+1, dim+1);

    R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
    R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
    R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
    R_inv_post(dim, dim) = temp_element;

    return R_inv_post.cwiseAbs().maxCoeff() < iter_colinear_threshold;
}
*/


bool calc_R_inverse_forward(
    const MatrixXd& R_pre,
    const VectorXd& r_temp_vec,
    double lower_right_corner,
    double iter_colinear_threshold,
    MatrixXd& R_post, 
    MatrixXd& R_inv_post,
    bool do_check,
    const ArrayXd& scaling_vector)
{   
    int dim = R_pre.rows();
    R_post.setZero(dim + 1, dim + 1);

    R_post.block(0, 0, dim, dim) = R_pre;
    R_post.block(0, dim, dim, 1) = r_temp_vec;
    R_post.block(dim, 0, 1, dim) = r_temp_vec.transpose();
    R_post(dim, dim) = lower_right_corner;

    LDLT<MatrixXd> ldlt_solver(R_post);
    ArrayXd ldlt_D = ldlt_solver.vectorD().array();
    R_inv_post = ldlt_solver.solve(MatrixXd::Identity(dim + 1, dim + 1));

    // sumstat_candidate.col(8) is X_norm_square, which has been set to 1 in read_sumstat() 
    // if if_keep_NA==false or if_LD_mode==false or if_gcta_COJO==false, then no changes to ldlt_D
    // the reason is that original gcta use unnormalized inner products for checking colinearity
    
    if (!do_check) return true;

    ldlt_D = ldlt_D / scaling_vector;

    return (ldlt_D.minCoeff() > 0) && (sqrt(ldlt_D.maxCoeff() / ldlt_D.minCoeff()) < 30) 
                && (R_inv_post.diagonal().maxCoeff() < iter_colinear_threshold);
}


// overload, remove one row and one column
// no need to check colinearity becuase it is guaranteed in forward step
void calc_R_inverse_backward(
    const MatrixXd& R_pre,
    const MatrixXd& R_inv_pre,
    int remove_index,
    MatrixXd& R_post,
    MatrixXd& R_inv_post)
{
    int dim = R_inv_pre.rows();
    if (remove_index < 0 || remove_index >= dim)
        throw std::invalid_argument("calc_block_inverse_fast: remove_index out of range.");
    /*
    if (if_fast_inv) {
        R_inv_post = R_inv_pre - R_inv_pre.col(remove_index) * R_inv_pre.row(remove_index) / R_inv_pre(remove_index, remove_index);
        remove_row(R_inv_post, remove_index);
        remove_column(R_inv_post, remove_index);
    }*/

    R_post = R_pre;
    remove_row(R_post, remove_index);
    remove_column(R_post, remove_index);
    R_inv_post = R_post.ldlt().solve(MatrixXd::Identity(dim - 1, dim - 1));
}


double median(const std::vector<double> &v)
{
    int size = v.size();
	if (size == 0) 
        throw std::invalid_argument("median: empty vector");

    std::vector<double> b(v);
    std::sort(b.begin(), b.end());
	return (size%2==1) ? b[(size-1)/2] : (b[size/2]+b[size/2-1])/2;
}


double median(const ArrayXd &eigen_vector)
{
    std::vector<double> v(eigen_vector.data(), eigen_vector.data() + eigen_vector.size());
    return median(v);
}
