#include "macojo.h"


long MACOJO::commonSNP_total_num;
map<string, int> MACOJO::commonSNP_index_map;
vector<string> MACOJO::A1_ref, MACOJO::A2_ref;
vector<int> MACOJO::SNP_pos_ref;


void Cohort::skim_cojo(string cojofile) 
{  
    ifstream Meta(cojofile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open the file [" + cojofile + "] to read");

    // SNP A1 A2 freq b se p N 
    string SNP_buf, str_buf;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line
    if (split_string(str_buf, vs_buf) < 7) LOGGER.e(0, "format error in the file [" + cojofile + "]");

    int SNP_num = 0;

    while (Meta) {
        Meta >> SNP_buf;
        if (Meta.eof()) break;
        Meta >> str_buf >> str_buf >> str_buf >> str_buf >> str_buf >> str_buf >> str_buf;
        SNP_cojo.push_back(SNP_buf);
        SNP_num++;
    }

    Meta.clear();
    Meta.close();

    sort(SNP_cojo.begin(), SNP_cojo.end());
    if (adjacent_find(SNP_cojo.begin(), SNP_cojo.end()) != SNP_cojo.end())
        LOGGER.e(0, "Duplicate SNP in [" + cojofile + "], please check");
    
    LOGGER << endl << SNP_num << " SNPs in sumstat file [" + cojofile + "]" << endl;
}


void Cohort::skim_PLINK(string PLINKfile) 
{   
    // Step 1: read .fam to get individual number 
    string famFile = PLINKfile+".fam";
    ifstream Fam(famFile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open FAM file [" + famFile + "] to read");

    string str_buf1, str_buf2;
    vector<string> indi_ids;
    indi_num = 0;

    while (Fam) {
        Fam >> str_buf1;
        if (Fam.eof()) break;
        Fam >> str_buf2;
        indi_ids.push_back(str_buf1+':'+str_buf2);
        Fam >> str_buf1 >> str_buf1 >> str_buf1 >> str_buf1;
        indi_num++;
    }

    Fam.clear();
    Fam.close();

    sort(indi_ids.begin(), indi_ids.end());
    if (adjacent_find(indi_ids.begin(), indi_ids.end()) != indi_ids.end())
        LOGGER.e(0, "Duplicate individual ID in FAM file [" + famFile + "], please check");
    
    // Step 2: read .bim to get SNP number
    string bimFile = PLINKfile+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "] to read");

    // Read bim file
    int SNP_num = 0;
    string SNP_buf, str_buf;

    while (Bim) {
        Bim >> str_buf;
        if (Bim.eof()) break;
        Bim >> SNP_buf >> str_buf >> str_buf >> str_buf >> str_buf;
        SNP_PLINK.push_back(SNP_buf);
        SNP_num++;
    }

    Bim.clear();
    Bim.close();
    
    sort(SNP_PLINK.begin(), SNP_PLINK.end());
    if (adjacent_find(SNP_PLINK.begin(), SNP_PLINK.end()) != SNP_PLINK.end())
        LOGGER.e(0, "Duplicate SNP in BIM file [" + bimFile + "], please check");

    LOGGER << indi_num << " individuals and " << SNP_num << " SNPs in PLINK file [" + PLINKfile + "]" << endl;
}


void Cohort::read_sumstat(string cojofile) 
{  
    LOGGER << endl << "Reading GWAS summary-level statistics from [" + cojofile + "] ..." << endl;
    ifstream Meta(cojofile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open the file [" + cojofile + "] to read");

    string SNP_buf, A1_buf, A2_buf, freq_buf, b_buf, se_buf, p_buf, N_buf, str_buf;
    double freq, b, se, p, N;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line
    if (split_string(str_buf, vs_buf) < 7) LOGGER.e(0, "format error in sumstat file [" + cojofile + "]");
    
    map<string, int>::iterator iter;
    int ref_index;

    sumstat.resize(MACOJO::commonSNP_total_num, 7);
    
    while (Meta) {
        Meta >> SNP_buf;
        if (Meta.eof()) break;
        Meta >> A1_buf >> A2_buf >> freq_buf >> b_buf >> se_buf >> p_buf >> N_buf;

        if ((iter = MACOJO::commonSNP_index_map.find(SNP_buf)) == MACOJO::commonSNP_index_map.end()) 
            continue;

        if (b_buf == "NA" || b_buf == "." || se_buf == "NA" || se_buf == "." || se_buf == "0" || \
            p_buf == "NA" || p_buf == "." || N_buf == "NA" || N_buf == "." || stod(N_buf) < 10) {
                LOGGER.w(1, "removed, invalid value in sumstat file [" + cojofile + "]", SNP_buf);
                MACOJO::commonSNP_index_map.erase(iter);
                continue;
            }

        to_upper(A1_buf);
        to_upper(A2_buf);
        freq = atof(freq_buf.c_str());
        b = atof(b_buf.c_str());
        se = atof(se_buf.c_str());
        p = atof(p_buf.c_str());
        N = atof(N_buf.c_str());

        ref_index = iter->second;

        if (MACOJO::A1_ref[ref_index] == A1_buf || MACOJO::A2_ref[ref_index] == A2_buf);
        else if (MACOJO::A1_ref[ref_index] == A2_buf && MACOJO::A2_ref[ref_index] == A1_buf) {freq = 1 - freq; b = -b;}
        else {
            LOGGER.w(1, "removed, A1 and A2 different from ref BIM file, please check", SNP_buf);
            MACOJO::commonSNP_index_map.erase(iter);
            continue;
        }

        // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
        sumstat(ref_index, 0) = b;
        sumstat(ref_index, 1) = se * se;
        sumstat(ref_index, 2) = p;
        sumstat(ref_index, 3) = freq;
        sumstat(ref_index, 4) = N;
    }

    Meta.close();
    LOGGER << "Finished reading sumstat file" << endl;
}


void Cohort::read_PLINK(string PLINKfile, bool is_ref_cohort) 
{   
    string bimFile = PLINKfile+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "] to read");
    LOGGER << endl << "Reading PLINK BIM file from [" + bimFile + "]..." << endl;

    string bedFile = PLINKfile+".bed";
    fstream Bed(bedFile.c_str(), ios::in | ios::binary);
    if (!Bed) LOGGER.e(0, "cannot open BED file [" + bedFile + "] to read");
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;

    // Read bim file
    int ibuf = 0;
    string SNP_buf, A1_buf = "0", A2_buf = "0", str_buf;

    // Read bed file
    int i, k;
    char ch[1];
    bitset<8> binary_buf;
    for (i=0; i<3; i++) {Bed.read(ch, 1);} // skip the first three bytes
    
    map<string, int>::iterator iter;
    long ref_index, X_start_index;

    bool SNP1, SNP2, swap;
    double SNP12, SNP_sum, SNP_square_sum, not_NA_indi_num, SNP_avg, SNP_std;
    
    X_A1.resize(MACOJO::commonSNP_total_num * indi_num);
    X_A2.resize(MACOJO::commonSNP_total_num * indi_num);
    X_avg.resize(MACOJO::commonSNP_total_num);
    X_std.resize(MACOJO::commonSNP_total_num);

    while (Bim) {
        Bim >> str_buf;
        if (Bim.eof()) break;
        Bim >> SNP_buf >> str_buf >> ibuf >> A1_buf >> A2_buf;

        if ((iter = MACOJO::commonSNP_index_map.find(SNP_buf)) == MACOJO::commonSNP_index_map.end()) {
            for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
            continue;
        }
        
        to_upper(A1_buf);
        to_upper(A2_buf);

        ref_index = iter->second;
        X_start_index = ref_index * indi_num;

        SNP_sum = 0; 
        SNP_square_sum = 0; 
        not_NA_indi_num = 0;
        swap = false;
        
        if (is_ref_cohort) {
            // read bim 1 as ref
            MACOJO::A1_ref[ref_index] = A1_buf;
            MACOJO::A2_ref[ref_index] = A2_buf;
            MACOJO::SNP_pos_ref[ref_index] = ibuf;
        } else {
            // check allele in bim 2
            if (MACOJO::A1_ref[ref_index] == A1_buf || MACOJO::A2_ref[ref_index] == A2_buf);
            else if (MACOJO::A1_ref[ref_index] == A2_buf && MACOJO::A2_ref[ref_index] == A1_buf)
                swap = true;
            else {
                LOGGER.w(1, "removed, A1 and A2 different between two BIM files, please check", SNP_buf);
                MACOJO::commonSNP_index_map.erase(iter);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }

            // check SNP position in bim 2
            if (ibuf != MACOJO::SNP_pos_ref[ref_index]) {
                LOGGER.w(1, "removed, SNP position different between two BIM files, please check", SNP_buf);
                MACOJO::commonSNP_index_map.erase(iter);
                for (i=0; i<indi_num; i+=4) Bed.read(ch, 1);
                continue;
            }
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

                if (swap && SNP1 == SNP2) {
                    SNP1 = !SNP1; SNP2 = !SNP2;
                }

                X_A1[X_start_index+i] = SNP1;
                X_A2[X_start_index+i] = SNP2;

                if (!SNP1 || SNP2) {
                    SNP12 = SNP1 + SNP2;
                    SNP_sum += SNP12;
                    SNP_square_sum += SNP12 * SNP12;
                    not_NA_indi_num++;
                }

                i++;
            }
        }
        
        if (not_NA_indi_num == 0) {
            LOGGER.w(1, "removed, all values are NA in bedfile", SNP_buf);
            MACOJO::commonSNP_index_map.erase(iter);
            continue;
        }

        SNP_avg = SNP_sum/not_NA_indi_num;
        SNP_std = sqrt((SNP_square_sum-SNP_avg*SNP_avg*not_NA_indi_num)/(indi_num-1));

        if (SNP_std < 1e-5) {
            LOGGER.w(1, "removed, identical genotypes for all individuals in bedfile", SNP_buf);
            MACOJO::commonSNP_index_map.erase(iter);
            continue;
        }
        
        X_avg[ref_index] = SNP_avg;
        X_std[ref_index] = SNP_std;
    }

    Bim.clear();
    Bim.close();
    Bed.clear();
    Bed.close();
    LOGGER << "Finished reading PLINK file" << endl;
}


void Cohort::get_vector_from_bed_matrix(int index, VectorXd &vec)
{   
    vec.setZero(indi_num);
    long X_start_index = long(index) * indi_num;

    #pragma omp parallel for
    for (int i = 0; i < indi_num; i++) {
        bool A1 = X_A1[X_start_index + i], A2 = X_A2[X_start_index + i];
        if (!A1 || A2) 
            vec(i) = (double(A1) + double(A2) - X_avg[index]) / X_std[index];
    }
}


void Cohort::calc_Vp(ArrayXXd &s) 
{   
    vector<string>().swap(SNP_cojo);
    vector<string>().swap(SNP_PLINK);

    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
    s.col(5) = s.col(3) * (1-s.col(3)) * 2;
    ArrayXd Vp_gcta_list = s.col(5) * s.col(4) * (s.col(1) + square(s.col(0))/(s.col(4)-1));
    Vp = median(Vp_gcta_list);

    s.col(4) = (Vp - s.col(5)*square(s.col(0))) / (s.col(5)*s.col(1)) + 1;
    s.col(6) = s.col(4) * s.col(5);

    Vp_gcta_list = s.col(6) * (s.col(1) + square(s.col(0))/(s.col(4)-1));
    Vp = median(Vp_gcta_list);
}


void MACOJO::read_files_two_cohorts(string cojoFile1, string PLINK1, string cojoFile2, string PLINK2) 
{   
    // scan all files to get rough common SNPs
    vector<string> commonSNP_c1, commonSNP_c2;

    // Cohort 1
    clock_t tStart = clock();
    c1.skim_cojo(cojoFile1);
    c1.skim_PLINK(PLINK1);
    set_intersection(c1.SNP_cojo.begin(), c1.SNP_cojo.end(), 
        c1.SNP_PLINK.begin(), c1.SNP_PLINK.end(), back_inserter(commonSNP_c1));
    LOGGER.i(0, "common SNPs included for Cohort 1", to_string(commonSNP_c1.size()));
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    
    // Cohort 2
    tStart = clock();
    c2.skim_cojo(cojoFile2);
    c2.skim_PLINK(PLINK2);
    set_intersection(c2.SNP_cojo.begin(), c2.SNP_cojo.end(), 
        c2.SNP_PLINK.begin(), c2.SNP_PLINK.end(), back_inserter(commonSNP_c2));
    LOGGER.i(0, "common SNPs included for Cohort 2", to_string(commonSNP_c2.size()));
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    
    // find common SNPs between two cohorts
    tStart = clock();
    auto first1 = commonSNP_c1.begin(), last1 = commonSNP_c1.end();
    auto first2 = commonSNP_c2.begin(), last2 = commonSNP_c2.end();
    int temp_index = 0;

    while (first1 != last1 && first2 != last2)
    {
        if (*first1 < *first2)  ++first1;
        else if (*first2 < *first1) ++first2;
        else {
            commonSNP_index_map.insert(commonSNP_index_map.end(), make_pair(*first1, temp_index));
            ++first1; ++first2; ++temp_index;
        }
    }
    
    vector<string>().swap(commonSNP_c1);
    vector<string>().swap(commonSNP_c2);
    
    A1_ref.resize(temp_index);
    A2_ref.resize(temp_index);
    SNP_pos_ref.resize(temp_index);
    commonSNP_total_num = temp_index;

    LOGGER << endl;
    LOGGER.i(0, "common SNPs included for two cohorts", to_string(temp_index));

    // set bim file 1 as reference
    c1.read_PLINK(PLINK1, true);
    LOGGER.i(0, "Remaining common SNPs", to_string(commonSNP_index_map.size()));
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    
    // read and adjust bim file 2 according to bim file 1
    tStart = clock();
    c2.read_PLINK(PLINK2, false);
    LOGGER.i(0, "Remaining common SNPs", to_string(commonSNP_index_map.size()));
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    
    // read sumstat files
    tStart = clock();
    c1.read_sumstat(cojoFile1);
    c2.read_sumstat(cojoFile2);

    // exclude rare SNP
    c1.sumstat_screened.setZero(commonSNP_index_map.size(), 7);
    c2.sumstat_screened.setZero(commonSNP_index_map.size(), 7);
    temp_index = 0;

    for (auto iter = commonSNP_index_map.begin(); iter != commonSNP_index_map.end(); ) {
        if (abs(c1.sumstat(iter->second, 3)-0.5) > 0.49 && abs(c2.sumstat(iter->second, 3)-0.5) > 0.49)
            iter = commonSNP_index_map.erase(iter);
        else {
            c1.sumstat_screened.row(temp_index) = c1.sumstat.row(iter->second);
            c2.sumstat_screened.row(temp_index) = c2.sumstat.row(iter->second);
            final_commonSNP.push_back(iter->first);
            final_commonSNP_index.push_back(iter->second);
            temp_index++;
            iter++;
        }
    }
    
    if (temp_index == 0)
        LOGGER.e(0, "Input data has no common SNPs");

    c1.sumstat_screened.conservativeResize(temp_index, 7);
    c1.calc_Vp(c1.sumstat_screened);
    c1.sumstat = c1.sumstat_screened;

    c2.sumstat_screened.conservativeResize(temp_index, 7);
    c2.calc_Vp(c2.sumstat_screened);
    c2.sumstat = c2.sumstat_screened;
    
    LOGGER << endl;
    LOGGER.i(0, "common SNPs at last", to_string(temp_index));
    LOGGER << "Vp: " << c1.Vp << " " << c2.Vp << endl;
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "--------------------------------" << endl << endl;
}


void MACOJO::read_files_one_cohort(string cojoFile, string PLINK) 
{   
    clock_t tStart = clock();
    c1.skim_cojo(cojoFile);
    c1.skim_PLINK(PLINK);

    auto first1 = c1.SNP_cojo.begin(), last1 = c1.SNP_cojo.end();
    auto first2 = c1.SNP_PLINK.begin(), last2 = c1.SNP_PLINK.end();
    int temp_index = 0;

    while (first1 != last1 && first2 != last2)
    {
        if (*first1 < *first2)  ++first1;
        else if (*first2 < *first1) ++first2;
        else {
            commonSNP_index_map.insert(commonSNP_index_map.end(), make_pair(*first1, temp_index));
            ++first1; ++first2; ++temp_index;
        }
    }

    LOGGER.i(0, "common SNPs included", to_string(temp_index));

    A1_ref.resize(temp_index);
    A2_ref.resize(temp_index);
    SNP_pos_ref.resize(temp_index);
    commonSNP_total_num = temp_index;

    // set bim file as reference
    c1.read_PLINK(PLINK, true);
    LOGGER.i(0, "Remaining common SNPs", to_string(commonSNP_index_map.size()));
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    
    // read sumstat files
    tStart = clock();
    c1.read_sumstat(cojoFile);

    // exclude rare SNP
    c1.sumstat_screened.setZero(commonSNP_index_map.size(), 7);
    temp_index = 0;

    for (auto iter = commonSNP_index_map.begin(); iter != commonSNP_index_map.end(); ) {
        if (abs(c1.sumstat(iter->second, 3)-0.5) > 0.49)
            iter = commonSNP_index_map.erase(iter);
        else {
            c1.sumstat_screened.row(temp_index) = c1.sumstat.row(iter->second);
            final_commonSNP.push_back(iter->first);
            final_commonSNP_index.push_back(iter->second);
            temp_index++;
            iter++; 
        }
    }

    if (temp_index == 0)
        LOGGER.e(0, "Input data has no common SNPs");

    c1.sumstat_screened.conservativeResize(temp_index, 7);
    c1.calc_Vp(c1.sumstat_screened);
    c1.sumstat = c1.sumstat_screened;

    LOGGER << endl;
    LOGGER.i(0, "common SNPs at last", to_string(temp_index));
    LOGGER << "Vp: " << c1.Vp << endl;
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "--------------------------------" << endl << endl;
}