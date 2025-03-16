#include "trans_ancestry_cojo.h"


void taCOJO::read_sumstat(string cojoFile, bool isFirst) 
{    
    LOGGER << "Reading GWAS summary-level statistics from [" + cojoFile + "] ..." << endl;
    ifstream Meta(cojoFile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open the file [" + cojoFile + "] to read");

    // SNP A1 A2 freq b se p N 
    double freq_buf = 0.0, b_buf = 0.0, se_buf = 0.0, p_buf = 0.0, N_buf = 0.0;
    string A1_buf, A2_buf, snp_buf, str_buf0, str_buf;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line
    if (StrFunc::split_string(str_buf, vs_buf) < 7) LOGGER.e(0, "format error in the input file [" + cojoFile + "]");

    int SNP_num = 0;

    while (Meta) {
        getline(Meta, str_buf0);
        stringstream iss(str_buf0);
        iss >> snp_buf >> A1_buf >> A2_buf;
        
        StrFunc::to_upper(A1_buf);
        StrFunc::to_upper(A2_buf);
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
        
        if (isFirst) {
            SNP_index_cohort1.insert(pair<string, int> (snp_buf, SNP_num));
            SNP_cohort1.push_back(snp_buf);
            A1_cohort1.push_back(A1_buf);
            A2_cohort1.push_back(A2_buf);
            freq1.push_back(freq_buf);
            b1.push_back(b_buf);
            se1.push_back(se_buf);
            p1.push_back(p_buf);
            N1.push_back(N_buf);
        } else {
            SNP_index_cohort2.insert(pair<string, int> (snp_buf, SNP_num));
            SNP_cohort2.push_back(snp_buf);
            A1_cohort2.push_back(A1_buf);
            A2_cohort2.push_back(A2_buf);
            freq2.push_back(freq_buf);
            b2.push_back(b_buf);
            se2.push_back(se_buf);
            p2.push_back(p_buf);
            N2.push_back(N_buf);
        }
        
        SNP_num++;
    }
    Meta.close();

    if ((isFirst && SNP_index_cohort1.size()<SNP_num) || (!isFirst && SNP_index_cohort2.size()<SNP_num))
        LOGGER.e(0, "Duplicate SNP in [" + cojoFile + "], please check");

    LOGGER << SNP_num << " SNPs included" << endl << endl;;
}


double taCOJO::read_PLINK(string PLINKfile, bool isFirst, vector<vector<double>> &gene, vector<string> &included_SNP_PLINK) 
{   
    string famFile = PLINKfile+".fam", bimFile = PLINKfile+".bim", bedFile = PLINKfile+".bed";

    // Step 1: read .fam to get individual number
    ifstream Fam(famFile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open FAM file [" + famFile + "].fam to read");
    LOGGER << "Reading PLINK FAM file from [" + famFile + "]..." << endl;

    string str_buf1, str_buf2, indi_id_buf;
    set<string> indi_ids;
    int indi_num = 0;

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
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "].bim to read");
    LOGGER << "Reading PLINK BIM file from [" + bimFile + "]..." << endl;

    fstream Bed(bedFile.c_str(), ios::in | ios::binary);
    if (!Bed) LOGGER.e(0, "cannot open BED file [" + bedFile + "].bed to read");
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;

    // Read bim file
    int ibuf = 0;
    string A1_buf = "0", A2_buf = "0";
    double dbuf = 0.0;

    // Read bed file
    int i, k;
    char ch[1];
    bitset<8> b;
    for (i=0; i<3; i++) {Bed.read(ch, 1);} // skip the first three bytes
    
    string SNP_buf;
    int SNP_buf_index;
    map<string, int>::iterator SNP_index_iter;

    bool SNP1, SNP2, swap;
    double SNP12, SNP_sum, SNP_square_sum, not_NA_indi_num, SNP_avg, SNP_std;
    vector<double> single_SNP_X_buf(indi_num);

    gene.resize(commonSNP.size());
    
    while (Bim) {
        Bim >> ibuf;
        if (Bim.eof()) break;
        Bim >> SNP_buf;
        Bim >> dbuf;
        Bim >> ibuf;
        Bim >> A1_buf;
        StrFunc::to_upper(A1_buf);
        Bim >> A2_buf;
        StrFunc::to_upper(A2_buf);

        if ((SNP_index_iter = commonSNP_index.find(SNP_buf)) == commonSNP_index.end()) {
            for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
            continue;
        }
        
        SNP_buf_index = SNP_index_iter->second;
        SNP_sum = 0;
        SNP_square_sum = 0;
        not_NA_indi_num = 0;
        swap = false;

        if (isFirst) {
            A1_ref[SNP_buf_index] = A1_buf;
            A2_ref[SNP_buf_index] = A2_buf;
            SNP_pos_ref[SNP_buf_index] = ibuf;

            if (abs(freq1[SNP_index_cohort1[SNP_buf]]-0.5) > 0.49 && abs(freq2[SNP_index_cohort2[SNP_buf]]-0.5) > 0.49) {
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;    
            }

            // check cohort 1 allele
            if (A1_buf == A2_cohort1[SNP_index_cohort1[SNP_buf]] && A2_buf == A1_cohort1[SNP_index_cohort1[SNP_buf]]) {
                freq1[SNP_index_cohort1[SNP_buf]] = 1-freq1[SNP_index_cohort1[SNP_buf]];
                b1[SNP_index_cohort1[SNP_buf]] = -b1[SNP_index_cohort1[SNP_buf]];
            } else if (A1_buf != A1_cohort1[SNP_index_cohort1[SNP_buf]] || A2_buf != A2_cohort1[SNP_index_cohort1[SNP_buf]]) {
                LOGGER.w(1, "SNP A1 and A2 different in sumstat file 1, please check", SNP_buf);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }

            // check cohort 2 allele
            if (A1_buf == A2_cohort2[SNP_index_cohort2[SNP_buf]] && A2_buf == A1_cohort2[SNP_index_cohort2[SNP_buf]]) {
                freq2[SNP_index_cohort2[SNP_buf]] = 1-freq2[SNP_index_cohort2[SNP_buf]];
                b2[SNP_index_cohort2[SNP_buf]] = -b2[SNP_index_cohort2[SNP_buf]];
            } else if (A1_buf != A1_cohort2[SNP_index_cohort2[SNP_buf]] || A2_buf != A2_cohort2[SNP_index_cohort2[SNP_buf]]) {
                LOGGER.w(1, "SNP A1 and A2 different in sumstat file 2, please check", SNP_buf);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }
        } else {
            // check SNP position in the second bim file
            if (ibuf != SNP_pos_ref[SNP_buf_index]) {
                LOGGER.w(1, "SNP position different in two BIM files, please check", SNP_buf);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }

            // check allele in the second bim file
            if (A2_buf == A1_ref[SNP_buf_index] && A1_buf == A2_ref[SNP_buf_index])
                swap = true;
            else if (A1_buf != A1_ref[SNP_buf_index] || A2_buf != A2_ref[SNP_buf_index]) {
                LOGGER.w(1, "SNP A1 and A2 different in two BIM files, please check", SNP_buf);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }
        }

        // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
        for (i = 0; i < indi_num;) {
            Bed.read(ch, 1);
            if (!Bed) LOGGER.e(0, "problem with the BED file ... has the FAM/BIM file been changed?");
            b = ch[0];
            k = 0;
            while (k < 7 && i < indi_num) {
                SNP2 = (!b[k++]);
                SNP1 = (!b[k++]);
                if (SNP1 && !SNP2) 
                    single_SNP_X_buf[i] = 10;
                else {
                    SNP12 = swap ? 2-SNP1-SNP2 : SNP1+SNP2;
                    single_SNP_X_buf[i] = SNP12;
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
        gene[SNP_buf_index].resize(indi_num);

        # pragma omp parallel for 
        for (i = 0; i < indi_num; i++) {
            gene[SNP_buf_index][i] = (single_SNP_X_buf[i] > 5) ? 0 : (single_SNP_X_buf[i]-SNP_avg)/SNP_std;
        }

        included_SNP_PLINK.push_back(SNP_buf);
    }

    Bim.clear();
    Bim.close();
    Bed.clear();
    Bed.close();

    vector<double>().swap(single_SNP_X_buf);

    sort(included_SNP_PLINK.begin(), included_SNP_PLINK.end());
    if (adjacent_find(included_SNP_PLINK.begin(), included_SNP_PLINK.end()) != included_SNP_PLINK.end())
        LOGGER.e(0, "Duplicate SNP in PLINK file [" + PLINKfile + "], please check");
    
    LOGGER << indi_num << " individuals and " << included_SNP_PLINK.size() << " SNPs included" << endl << endl;
    return indi_num;
}


void taCOJO::generate_sumstat_and_X() 
{   
    // initialize X matrices
    // row: individual, column: SNP
    
    X1.resize(indi_num1, commonSNP_num);
    X2.resize(indi_num2, commonSNP_num);

    # pragma omp parallel for 
    for (int i = 0; i < commonSNP_num; i++) {
        int index = commonSNP_index[commonSNP[i]];
        X1.col(i) = Map<ArrayXd>(gene_cohort1[index].data(), indi_num1);
        X2.col(i) = Map<ArrayXd>(gene_cohort2[index].data(), indi_num2);
    }
    
    vector<vector<double>>().swap(gene_cohort1);
    vector<vector<double>>().swap(gene_cohort2);
    
    // initialize sumstat matrices
    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 

    sumstat1.resize(commonSNP_num, 7);
    sumstat2.resize(commonSNP_num, 7);

    # pragma omp parallel for 
    for (int i = 0; i < commonSNP_num; i++) {
        int index = commonSNP_index[commonSNP[i]];

        int cohort1_index = SNP_index_cohort1[commonSNP[i]];
        sumstat1(i, 0) = b1[cohort1_index];
        sumstat1(i, 1) = se1[cohort1_index] * se1[cohort1_index];
        sumstat1(i, 2) = p1[cohort1_index];
        sumstat1(i, 3) = freq1[cohort1_index];
        sumstat1(i, 4) = N1[cohort1_index];

        int cohort2_index = SNP_index_cohort2[commonSNP[i]];
        sumstat2(i, 0) = b2[cohort2_index];
        sumstat2(i, 1) = se2[cohort2_index] * se2[cohort2_index];
        sumstat2(i, 2) = p2[cohort2_index];
        sumstat2(i, 3) = freq2[cohort2_index];
        sumstat2(i, 4) = N2[cohort2_index];
    }

    map<string, int>().swap(SNP_index_cohort1);
    map<string, int>().swap(SNP_index_cohort2);
    
    vector<double>().swap(freq1);
    vector<double>().swap(b1);
    vector<double>().swap(se1);
    vector<double>().swap(p1);
    vector<double>().swap(N1);

    vector<double>().swap(freq2);
    vector<double>().swap(b2);
    vector<double>().swap(se2);
    vector<double>().swap(p2);
    vector<double>().swap(N2);
}


double taCOJO::median(const ArrayXd &eigen_vector) 
{
    if (eigen_vector.size()==1) 
        return eigen_vector(0);

    int size = eigen_vector.size();
    vector<double> b(eigen_vector.data(), eigen_vector.data() + size);
    double b_median; 

    stable_sort(b.begin(), b.end());
    if (size%2==1)
        b_median = b[(size-1)/2];
    else 
        b_median = (b[size/2]+b[size/2-1])/2;

    vector<double>().swap(b);
    return b_median;
}


double taCOJO::calc_Vp(ArrayXXd &sumstat) 
{
    sumstat.col(5) = sumstat.col(3) * (1-sumstat.col(3)) * 2;
    ArrayXd Vp_gcta_list = sumstat.col(5) * sumstat.col(4) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    double Vp = median(Vp_gcta_list);

    sumstat.col(4) = (Vp - sumstat.col(5)*square(sumstat.col(0))) / (sumstat.col(5)*sumstat.col(1)) + 1;
    sumstat.col(6) = sumstat.col(4) * sumstat.col(5);

    Vp_gcta_list = sumstat.col(6) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    Vp = median(Vp_gcta_list);    
    return Vp;
}


void taCOJO::read_files(string cojoFile1, string cojoFile2, string PLINK1, string PLINK2) 
{   
    // read cojofiles and get rough common SNPs
    read_sumstat(cojoFile1, true);
    read_sumstat(cojoFile2, false);

    sort(SNP_cohort1.begin(), SNP_cohort1.end());
    sort(SNP_cohort2.begin(), SNP_cohort2.end());
    set_intersection(SNP_cohort1.begin(), SNP_cohort1.end(), 
        SNP_cohort2.begin(), SNP_cohort2.end(), back_inserter(commonSNP));
        
    vector<string>().swap(SNP_cohort1); 
    vector<string>().swap(SNP_cohort2); 

    int i = 0;
    for (auto iter = commonSNP.begin(); iter != commonSNP.end(); iter++, i++)
        commonSNP_index.insert(commonSNP_index.end(), pair<string, int> (*iter, i));

    // set bim file 1 as reference, compare and get common SNPs across 4 files
    A1_ref.resize(commonSNP.size());
    A2_ref.resize(commonSNP.size());
    SNP_pos_ref.resize(commonSNP.size());

    indi_num1 = read_PLINK(PLINK1, true, gene_cohort1, included_SNP_PLINK_cohort1);
    indi_num2 = read_PLINK(PLINK2, false, gene_cohort2, included_SNP_PLINK_cohort2);
    
    vector<string>().swap(A1_cohort1);
    vector<string>().swap(A2_cohort1);
    vector<string>().swap(A1_cohort2);
    vector<string>().swap(A2_cohort2);

    if (commonSNP.size() != included_SNP_PLINK_cohort1.size() || 
        commonSNP.size() != included_SNP_PLINK_cohort2.size()) {

        vector<string>().swap(commonSNP);
        set_intersection(included_SNP_PLINK_cohort1.begin(), included_SNP_PLINK_cohort1.end(), 
            included_SNP_PLINK_cohort2.begin(), included_SNP_PLINK_cohort2.end(), back_inserter(commonSNP));
        
        vector<int> SNP_pos_ref_temp;

        for (auto iter = commonSNP.begin(); iter != commonSNP.end(); iter++)
            SNP_pos_ref_temp.push_back(SNP_pos_ref[commonSNP_index[*iter]]);
        
        SNP_pos_ref_temp.swap(SNP_pos_ref);
        vector<int>().swap(SNP_pos_ref_temp);        
    }

    vector<string>().swap(included_SNP_PLINK_cohort1); 
    vector<string>().swap(included_SNP_PLINK_cohort2); 

    commonSNP_num = commonSNP.size();
    if (commonSNP_num==0)
        LOGGER.e(0, "Input data has no common SNPs.");

    LOGGER.i(0, "common SNPs included in Genotype data", to_string(commonSNP_num));
    LOGGER.i(0, "individuals in Cohort 1", to_string(indi_num1));
    LOGGER.i(0, "individuals in Cohort 2", to_string(indi_num2));

    // initialize sumstat and X matrices, calculate Vp
    generate_sumstat_and_X();
    Vp1 = calc_Vp(sumstat1);
    Vp2 = calc_Vp(sumstat2);

    cout << "Vp: " << Vp1 << " " << Vp2 << endl << endl;
}
