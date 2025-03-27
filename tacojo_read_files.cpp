#include "trans_ancestry_cojo.h"


vector<string> Cohort::commonSNP;
map<string, int> Cohort::sumstat_commonSNP_index;
vector<string> Cohort::A1_ref, Cohort::A2_ref;
vector<int> Cohort::SNP_pos_ref;


void Cohort::read_sumstat(string cojofile) 
{  
    LOGGER << "Reading GWAS summary-level statistics from [" + cojofile + "] ..." << endl;
    ifstream Meta(cojofile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open the file [" + cojofile + "] to read");

    // SNP A1 A2 freq b se p N 
    double freq_buf = 0.0, b_buf = 0.0, se_buf = 0.0, p_buf = 0.0, N_buf = 0.0;
    string A1_buf, A2_buf, snp_buf, str_buf0, str_buf;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line
    if (split_string(str_buf, vs_buf) < 7) LOGGER.e(0, "format error in the input file [" + cojofile + "]");

    int SNP_num = 0;

    while (Meta) {
        getline(Meta, str_buf0);
        stringstream iss(str_buf0);
        iss >> snp_buf >> A1_buf >> A2_buf;
        
        to_upper(A1_buf);
        to_upper(A2_buf);
        iss >> str_buf;
        freq_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        b_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == "." || str_buf == "0") continue;
        se_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        p_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        N_buf = atof(str_buf.c_str());
        if (N_buf < 10) LOGGER.e(0, "invalid sample size in line:\n\"" + str_buf0 + "\"");
        if (Meta.eof()) break;
        
        SNP_index.insert(pair<string, int> (snp_buf, SNP_num));
        SNP.push_back(snp_buf);
        A1.push_back(A1_buf);
        A2.push_back(A2_buf);
        freq.push_back(freq_buf);
        b.push_back(b_buf);
        se.push_back(se_buf);
        p.push_back(p_buf);
        N.push_back(N_buf);
        
        SNP_num++;
    }
    Meta.close();

    if (SNP_index.size() < SNP_num)
        LOGGER.e(0, "Duplicate SNP in [" + cojofile + "], please check");
    
    sort(SNP.begin(), SNP.end());
    LOGGER << SNP_num << " SNPs included" << endl << endl;
}


void Cohort::read_PLINK(string PLINKfile, bool is_ref_cohort) 
{   
    // Step 1: read .fam to get individual number 
    string famFile = PLINKfile+".fam";
    ifstream Fam(famFile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open FAM file [" + famFile + "].fam to read");
    LOGGER << "Reading PLINK FAM file from [" + famFile + "]..." << endl;

    string str_buf1, str_buf2, indi_id_buf;
    set<string> indi_ids;
    indi_num = 0;

    while (Fam) {
        Fam >> str_buf1;
        if (Fam.eof()) break;
        Fam >> str_buf2;
        indi_ids.insert(str_buf1+':'+str_buf2);
        Fam >> str_buf1;
        Fam >> str_buf1;
        Fam >> str_buf1;
        Fam >> str_buf1;
        indi_num++;
    }
    Fam.clear();
    Fam.close();

    if (indi_num != indi_ids.size())
        LOGGER.e(0, "Duplicate individual ID in [" + famFile + "], please check");
    
    set<string>().swap(indi_ids);

    // Step 2: read .bim and .bed to get X matrix
    string bimFile = PLINKfile+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "].bim to read");
    LOGGER << "Reading PLINK BIM file from [" + bimFile + "]..." << endl;

    string bedFile = PLINKfile+".bed";
    fstream Bed(bedFile.c_str(), ios::in | ios::binary);
    if (!Bed) LOGGER.e(0, "cannot open BED file [" + bedFile + "].bed to read");
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;

    // Read bim file
    int ibuf = 0;
    string SNP_buf, A1_buf = "0", A2_buf = "0";
    double dbuf = 0.0;

    // Read bed file
    int i, k;
    char ch[1];
    bitset<8> binary_buf;
    for (i=0; i<3; i++) {Bed.read(ch, 1);} // skip the first three bytes
    
    map<string, int>::iterator iter;
    int ref_index, cohort_index;
    
    string A1_cohort, A2_cohort, A1_ref_buf, A2_ref_buf;

    bool SNP1, SNP2, swap;
    double SNP12, SNP_sum, SNP_square_sum, not_NA_indi_num, SNP_avg, SNP_std;
    vector<double> single_SNP_buf(indi_num);

    genotype.resize(sumstat_commonSNP_index.size());

    while (Bim) {
        Bim >> ibuf;
        if (Bim.eof()) break;
        Bim >> SNP_buf;
        Bim >> dbuf;
        Bim >> ibuf;
        Bim >> A1_buf;
        Bim >> A2_buf;
        to_upper(A1_buf);
        to_upper(A2_buf);

        if ((iter = sumstat_commonSNP_index.find(SNP_buf)) == sumstat_commonSNP_index.end()) {
            for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
            continue;
        }
        
        ref_index = iter->second;
        cohort_index = SNP_index[SNP_buf];

        SNP_sum = 0;
        SNP_square_sum = 0;
        not_NA_indi_num = 0;
        swap = false;

        A1_cohort = A1[cohort_index];
        A2_cohort = A2[cohort_index];
        
        if (is_ref_cohort) {
            A1_ref_buf = A1_buf;
            A2_ref_buf = A2_buf;

            // read bim 1 as ref
            A1_ref[ref_index] = A1_buf;
            A2_ref[ref_index] = A2_buf;
            SNP_pos_ref[ref_index] = ibuf;
        } else {
            A1_ref_buf = A1_ref[ref_index];
            A2_ref_buf = A2_ref[ref_index];

            // check allele in bim 2
            if (A1_ref_buf == A1_buf || A2_ref_buf == A2_buf);
            else if (A1_ref_buf == A2_buf && A2_ref_buf == A1_buf)
                swap = true;
            else {
                LOGGER.w(1, "A1 and A2 different between two BIM files, please check", SNP_buf);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }

            // check SNP position in bim 2
            if (ibuf != SNP_pos_ref[ref_index]) {
                LOGGER.w(1, "SNP position different between two BIM files, please check", SNP_buf);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }
        }
        
        // check allele in sumstat file
        if (A1_ref_buf == A1_cohort || A2_ref_buf == A2_cohort);
        else if (A1_ref_buf == A2_cohort && A2_ref_buf == A1_cohort) {
            freq[cohort_index] = 1-freq[cohort_index];
            b[cohort_index] = -b[cohort_index];
        } else {
            LOGGER.w(1, "A1 and A2 different between [" + bimFile + "] and sumstat file, please check", SNP_buf);
            for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
            continue;
        }

        // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
        for (i = 0; i < indi_num;) {
            Bed.read(ch, 1);
            if (!Bed) LOGGER.e(0, "problem with the BED file ... has the FAM/BIM file been changed?");
            binary_buf = ch[0];
            k = 0;
            while (k < 7 && i < indi_num) {
                SNP2 = (!binary_buf[k++]);
                SNP1 = (!binary_buf[k++]);
                if (SNP1 && !SNP2) 
                    single_SNP_buf[i] = 10;
                else {
                    SNP12 = swap ? 2-SNP1-SNP2 : SNP1+SNP2;
                    single_SNP_buf[i] = SNP12;
                    SNP_sum += SNP12;
                    SNP_square_sum += SNP12*SNP12;
                    not_NA_indi_num++;
                }
                i++;
            }
        }
        
        if (not_NA_indi_num == 0) {
            LOGGER.w(1, "all values are NA in bedfile " + bedFile, SNP_buf);
            continue;
        }

        SNP_avg = SNP_sum/not_NA_indi_num;
        SNP_std = sqrt((SNP_square_sum-SNP_avg*SNP_avg*not_NA_indi_num)/(indi_num-1));

        if (SNP_std < 1e-5) {
            LOGGER.w(1, "all genotypes are identical in bedfile " + bedFile, SNP_buf);
            continue;
        }
        
        // fill NA with 0
        genotype[ref_index].resize(indi_num);

        # pragma omp parallel for 
        for (i = 0; i < indi_num; i++) {
            genotype[ref_index][i] = (single_SNP_buf[i] > 5) ? 0 : (single_SNP_buf[i]-SNP_avg)/SNP_std;
        }

        included_SNP_PLINK.push_back(SNP_buf);
    }
    Bim.clear();
    Bim.close();
    Bed.clear();
    Bed.close();

    vector<double>().swap(single_SNP_buf);

    sort(included_SNP_PLINK.begin(), included_SNP_PLINK.end());
    if (adjacent_find(included_SNP_PLINK.begin(), included_SNP_PLINK.end()) != included_SNP_PLINK.end())
        LOGGER.e(0, "Duplicate SNP in PLINK file [" + PLINKfile + "], please check");
    
    LOGGER << indi_num << " individuals and " << included_SNP_PLINK.size() << " SNPs included" << endl << endl;
}


void Cohort::generate_sumstat_and_X() 
{   
    int commonSNP_num = commonSNP.size();
    
    X.resize(indi_num, commonSNP_num);
    sumstat.resize(commonSNP_num, 7);
    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 

    # pragma omp parallel for 
    for (int i = 0; i < commonSNP_num; i++) {
        int ref_index = sumstat_commonSNP_index[commonSNP[i]];
        int cohort_index = SNP_index[commonSNP[i]];

        X.col(i) = Map<ArrayXd>(genotype[ref_index].data(), indi_num);
        sumstat(i, 0) = b[cohort_index];
        sumstat(i, 1) = se[cohort_index]*se[cohort_index];
        sumstat(i, 2) = p[cohort_index];
        sumstat(i, 3) = freq[cohort_index];
        sumstat(i, 4) = N[cohort_index];
    }

    // release memory
    map<string, int>().swap(SNP_index);
    vector<string>().swap(A1);
    vector<string>().swap(A2);
    vector<double>().swap(freq);
    vector<double>().swap(b);
    vector<double>().swap(se);
    vector<double>().swap(p);
    vector<double>().swap(N);
    vector<vector<double>> genotype;
    vector<string>().swap(SNP);
    vector<string>().swap(included_SNP_PLINK);

    // calculate Vp
    sumstat.col(5) = sumstat.col(3) * (1-sumstat.col(3)) * 2;
    ArrayXd Vp_gcta_list = sumstat.col(5) * sumstat.col(4) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    Vp = median(Vp_gcta_list);

    sumstat.col(4) = (Vp - sumstat.col(5)*square(sumstat.col(0))) / (sumstat.col(5)*sumstat.col(1)) + 1;
    sumstat.col(6) = sumstat.col(4) * sumstat.col(5);

    Vp_gcta_list = sumstat.col(6) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    Vp = median(Vp_gcta_list);
}


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


void Cohort::calc_R_inv_fast() {
    double temp_element = 1 / (1 - r_temp_vec.transpose() * R_inv_pre * r_temp_vec);
    VectorXd temp_vector = R_inv_pre * r_temp_vec;

    int dim = R_inv_pre.rows();
    R_inv_post.resize(dim+1, dim+1);

    R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
    R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
    R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
    R_inv_post(dim, dim) = temp_element;
}


double median(const ArrayXd &eigen_vector)
{
    int size = eigen_vector.size();
    vector<double> b(eigen_vector.data(), eigen_vector.data() + size);
    double b_median; 

    sort(b.begin(), b.end());
    if (size%2==1)
        b_median = b[(size-1)/2];
    else 
        b_median = (b[size/2]+b[size/2-1])/2;

    vector<double>().swap(b);
    return b_median;
}


void to_upper(string &str)
{
	int i=0;
	for(i=0; i<str.size(); i++){
		if(str[i]>='a' && str[i]<='z') str[i]+='A'-'a';
	}
}


int split_string(const string &str, vector<string> &vec_str, string separator)
{
	if(str.empty()) return 0;
	vec_str.clear();

	int i=0;
	bool look=false;
	string str_buf;
	string symbol_pool="`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
	string::size_type pos;

	for(i=0; i<separator.size(); i++){
		pos=symbol_pool.find(separator[i]);
		if( pos!=string::npos ) symbol_pool.erase(symbol_pool.begin()+pos);
	}

	for(i=0; i<str.size(); i++){
		if( symbol_pool.find(str[i])!=string::npos ){
			if(!look) look=true;
			str_buf += str[i];
		}
		else{
			if(look){
				look=false;
				vec_str.push_back(str_buf);
				str_buf.erase(str_buf.begin(), str_buf.end());
			}
		}
	}
	if(look) vec_str.push_back(str_buf);

	return vec_str.size();
}


void append_row(ArrayXXd &matrix, const ArrayXXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void append_row(MatrixXd &matrix, const MatrixXd &vector)
{   
    int numRows = matrix.rows();
    matrix.conservativeResize(numRows+1, NoChange);
    matrix.row(numRows) = vector;
}


void append_column(MatrixXd &matrix, const MatrixXd &vector)
{
    int numCols = matrix.cols();
    matrix.conservativeResize(NoChange, numCols+1);
    matrix.col(numCols) = vector;
}


void remove_row(ArrayXXd &matrix, int index)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void remove_row(MatrixXd &matrix, int index)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.middleRows(index, numRows-index) = matrix.bottomRows(numRows-index).eval();

    matrix.conservativeResize(numRows, NoChange);
}


void remove_column(MatrixXd &matrix, int index)
{   
    // -1 indicates the last column
    int numRows = matrix.rows(), numCols = matrix.cols()-1;

    if (index != -1 && index < numCols)
        matrix.middleCols(index, numCols-index) = matrix.rightCols(numCols-index).eval();

    matrix.conservativeResize(NoChange, numCols);
}
