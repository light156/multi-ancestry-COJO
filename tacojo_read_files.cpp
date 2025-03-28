#include "trans_ancestry_cojo.h"


vector<string> Cohort::commonSNP;
map<string, int> Cohort::sumstat_commonSNP_index;
vector<string> Cohort::A1_ref, Cohort::A2_ref;
vector<int> Cohort::SNP_pos_ref;


void Cohort::read_sumstat(string cojofile) 
{  
    LOGGER << endl << "Reading GWAS summary-level statistics from [" + cojofile + "] ..." << endl;
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
    LOGGER << SNP_num << " SNPs included" << endl;
}


void Cohort::read_PLINK(string PLINKfile, bool is_ref_cohort) 
{   
    // Step 1: read .fam to get individual number 
    string famFile = PLINKfile+".fam";
    ifstream Fam(famFile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open FAM file [" + famFile + "] to read");
    LOGGER << endl << "Reading PLINK FAM file from [" + famFile + "]..." << endl;

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
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "] to read");
    LOGGER << "Reading PLINK BIM file from [" + bimFile + "]..." << endl;

    string bedFile = PLINKfile+".bed";
    fstream Bed(bedFile.c_str(), ios::in | ios::binary);
    if (!Bed) LOGGER.e(0, "cannot open BED file [" + bedFile + "] to read");
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
    vector<string>().swap(A1);
    vector<string>().swap(A2);

    sort(included_SNP_PLINK.begin(), included_SNP_PLINK.end());
    if (adjacent_find(included_SNP_PLINK.begin(), included_SNP_PLINK.end()) != included_SNP_PLINK.end())
        LOGGER.e(0, "Duplicate SNP in PLINK file [" + PLINKfile + "], please check");
    
    LOGGER << indi_num << " individuals and " << included_SNP_PLINK.size() << " SNPs included" << endl;
}


void Cohort::finalize_ref_info() 
{
    vector<string> A1_ref_clean, A2_ref_clean;
    vector<int> SNP_pos_ref_clean;
    int temp_index;

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
}


void Cohort::generate_sumstat_and_X() 
{   
    int commonSNP_num = commonSNP.size();
    
    sumstat.resize(commonSNP_num, 7);
    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 

    # pragma omp parallel for 
    for (int i = 0; i < commonSNP_num; i++) {
        int cohort_index = SNP_index[commonSNP[i]];
        sumstat(i, 0) = b[cohort_index];
        sumstat(i, 1) = se[cohort_index]*se[cohort_index];
        sumstat(i, 2) = p[cohort_index];
        sumstat(i, 3) = freq[cohort_index];
        sumstat(i, 4) = N[cohort_index];
    }

    // release memory
    map<string, int>().swap(SNP_index);
    vector<double>().swap(freq);
    vector<double>().swap(b);
    vector<double>().swap(se);
    vector<double>().swap(p);
    vector<double>().swap(N);

    X.resize(indi_num, commonSNP_num);

    # pragma omp parallel for 
    for (int i = 0; i < commonSNP_num; i++) {
        int ref_index = sumstat_commonSNP_index[commonSNP[i]];
        X.col(i) = Map<ArrayXd>(genotype[ref_index].data(), indi_num);
    }

    // release memory
    vector<vector<double>> genotype;

    // calculate Vp
    sumstat.col(5) = sumstat.col(3) * (1-sumstat.col(3)) * 2;
    ArrayXd Vp_gcta_list = sumstat.col(5) * sumstat.col(4) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    Vp = median(Vp_gcta_list);

    sumstat.col(4) = (Vp - sumstat.col(5)*square(sumstat.col(0))) / (sumstat.col(5)*sumstat.col(1)) + 1;
    sumstat.col(6) = sumstat.col(4) * sumstat.col(5);

    Vp_gcta_list = sumstat.col(6) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    Vp = median(Vp_gcta_list);
}


void TransAncestryCOJO::read_files_two_cohorts(string cojoFile1, string PLINK1, string cojoFile2, string PLINK2) 
{   
    clock_t tStart = clock();

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

    vector<string>().swap(c1.SNP);
    vector<string>().swap(c2.SNP);
    vector<string>().swap(Cohort::commonSNP);

    // set bim file 1 as reference, compare and get common SNPs across 4 files
    Cohort::A1_ref.resize(temp_index);
    Cohort::A2_ref.resize(temp_index);
    Cohort::SNP_pos_ref.resize(temp_index);

    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    tStart = clock();

    c1.read_PLINK(PLINK1, true);

    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    tStart = clock();

    c2.read_PLINK(PLINK2, false);
    
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    tStart = clock();

    set_intersection(c1.included_SNP_PLINK.begin(), c1.included_SNP_PLINK.end(), 
        c2.included_SNP_PLINK.begin(), c2.included_SNP_PLINK.end(), back_inserter(Cohort::commonSNP));
    
    vector<string>().swap(c1.included_SNP_PLINK);
    vector<string>().swap(c2.included_SNP_PLINK);

    if (Cohort::commonSNP.size() == 0)
        LOGGER.e(0, "Input data has no common SNPs.");

    LOGGER.i(0, "common SNPs included for two cohorts", to_string(Cohort::commonSNP.size()));
    LOGGER.i(0, "individuals in Cohort 1", to_string(c1.indi_num));
    LOGGER.i(0, "individuals in Cohort 2", to_string(c2.indi_num));

    // initialize sumstat and X matrices, calculate Vp
    c1.generate_sumstat_and_X();
    c2.generate_sumstat_and_X();
    LOGGER << endl << "Vp: " << c1.Vp << " " << c2.Vp << endl;

    // clean ref A1, A2 and SNP position lists
    Cohort::finalize_ref_info();
    
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "--------------------------------" << endl << endl;
}


void TransAncestryCOJO::read_files_one_cohort(string cojoFile1, string PLINK) 
{   
    clock_t tStart = clock();

    // read cojofiles and get rough common SNPs    
    c1.read_sumstat(cojoFile1);
    
    // exclude rare SNP
    int temp_index = 0;
    for (auto iter = c1.SNP.begin(); iter != c1.SNP.end(); iter++) {
        if (abs(c1.freq[c1.SNP_index[*iter]]-0.5) <= 0.49) {
            Cohort::sumstat_commonSNP_index.insert(Cohort::sumstat_commonSNP_index.end(), pair<string, int> (*iter, temp_index));
            temp_index++;
        }
    }

    vector<string>().swap(c1.SNP);

    // set bim file 1 as reference, compare and get common SNPs across 4 files
    Cohort::A1_ref.resize(temp_index);
    Cohort::A2_ref.resize(temp_index);
    Cohort::SNP_pos_ref.resize(temp_index);

    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    tStart = clock();

    c1.read_PLINK(PLINK, true);

    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    tStart = clock();
    
    c1.included_SNP_PLINK.swap(Cohort::commonSNP);

    if (Cohort::commonSNP.size() == 0)
        LOGGER.e(0, "Input data has no SNPs.");

    // initialize sumstat and X matrices, calculate Vp
    c1.generate_sumstat_and_X();
    LOGGER << endl << "Vp: " << c1.Vp << endl;

    // clean ref A1, A2 and SNP position lists
    Cohort::finalize_ref_info();

    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "--------------------------------" << endl << endl;
}
