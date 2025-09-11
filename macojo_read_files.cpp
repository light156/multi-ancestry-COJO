#include "macojo.h"


long MACOJO::commonSNP_total_num;
map<string, int> MACOJO::commonSNP_index_map;
vector<string> MACOJO::A1_ref, MACOJO::A2_ref;
vector<int> MACOJO::SNP_pos_ref;
vector<int> MACOJO::final_commonSNP_index;


void MACOJO::read_SNP_only(string filename, vector<string> &SNP_list, bool table_head) 
{
    ifstream Meta(filename.c_str());
    if (!Meta) LOGGER.e(0, "cannot open file [" + filename + "] to read");

    string str_buf, SNP_buf, temp;
    int SNP_num = 0;

    if (table_head) getline(Meta, str_buf); // the header line

    while (Meta) {
        getline(Meta, str_buf);
        if (Meta.eof()) break;
        if (str_buf.empty()) continue; // skip empty lines

        istringstream iss(str_buf);
        iss >> SNP_buf;
        SNP_list.push_back(SNP_buf);
        SNP_num++;
    }

    Meta.clear();
    Meta.close();

    sort(SNP_list.begin(), SNP_list.end());
    if (adjacent_find(SNP_list.begin(), SNP_list.end()) != SNP_list.end())
        LOGGER.e(0, "Duplicate SNP in [" + filename + "], please check", *adjacent_find(SNP_list.begin(), SNP_list.end()));
    
    LOGGER << endl << SNP_num << " SNPs in file [" + filename + "]" << endl;
}


void Cohort::skim_fam(string famFile) 
{   
    // read .fam to get individual number 
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

    LOGGER << indi_num << " individuals in FAM file [" + famFile + "]" << endl;
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
            // LOGGER.w(1, "removed, identical genotypes for all individuals in bedfile", SNP_buf);
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


void Cohort::calc_Vp() 
{   
    ArrayXXd &s = sumstat_screened;

    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
    s.col(5) = s.col(3) * (1-s.col(3)) * 2;
    ArrayXd Vp_gcta_list = s.col(5) * s.col(4) * (s.col(1) + square(s.col(0))/(s.col(4)-1));
    Vp = median(Vp_gcta_list);

    s.col(4) = (Vp - s.col(5)*square(s.col(0))) / (s.col(5)*s.col(1)) + 1;
    s.col(6) = s.col(4) * s.col(5);

    Vp_gcta_list = s.col(6) * (s.col(1) + square(s.col(0))/(s.col(4)-1));
    Vp = median(Vp_gcta_list);
}


void MACOJO::read_cojo_PLINK_files(char** filenames, int cohort_num) 
{   
    // Step 1: get common SNPs across all cohorts
    vector<string> commonSNP_all_cohorts;
    clock_t tStart;
    int temp_index = 0;

    for (int n=0; n<cohort_num; n++) {

        tStart = clock();

        auto &c = cohorts[n];
        vector<string> SNP_cojo, SNP_PLINK, SNP_common;

        string cojo_file = filenames[n*2], PLINK_file = filenames[n*2+1];

        read_SNP_only(cojo_file, SNP_cojo, true);
        read_SNP_only(PLINK_file+".bim", SNP_PLINK, false);
        c.skim_fam(PLINK_file+".fam");

        set_intersection(SNP_cojo.begin(), SNP_cojo.end(), 
            SNP_PLINK.begin(), SNP_PLINK.end(), back_inserter(SNP_common));

        if (n==0)
            commonSNP_all_cohorts = SNP_common;
        else {
            vector<string> temp;
            set_intersection(commonSNP_all_cohorts.begin(), commonSNP_all_cohorts.end(), 
                SNP_common.begin(), SNP_common.end(), back_inserter(temp));
            commonSNP_all_cohorts.swap(temp);
        }

        LOGGER << commonSNP_all_cohorts.size() << " common SNPs after including Cohort " << n+1 << endl;
        LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    } 

    if (all_SNP.size() > 0) {
        vector<string> temp;
        set_intersection(all_SNP.begin(), all_SNP.end(), 
            commonSNP_all_cohorts.begin(), commonSNP_all_cohorts.end(), back_inserter(temp));
        commonSNP_all_cohorts.swap(temp);
        LOGGER << "After user-given SNP list, common SNPs: " << commonSNP_all_cohorts.size() << endl;
    }
                       
    tStart = clock();

    for (const auto &snp : commonSNP_all_cohorts) {
        commonSNP_index_map.insert(make_pair(snp, temp_index));
        temp_index++;
    }
    
    A1_ref.resize(temp_index);
    A2_ref.resize(temp_index);
    SNP_pos_ref.resize(temp_index);
    commonSNP_total_num = temp_index;

    // Step 2: read bim and bed files for all cohorts
    for (int n=0; n<cohort_num; n++) {
        // set bim file 1 as reference
        tStart = clock();
        cohorts[n].read_PLINK(filenames[n*2+1], n==0);
        LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    }
    
    // Step 3: read sumstat files for all cohorts
    tStart = clock();

    for (int n=0; n<cohort_num; n++) 
        cohorts[n].read_sumstat(filenames[n*2]);

    // exclude rare SNP  
    for (auto &c : cohorts) 
        c.sumstat_screened.setZero(commonSNP_index_map.size(), 7);

    temp_index = 0;

    for (auto iter = commonSNP_index_map.begin(); iter != commonSNP_index_map.end(); ) {
        bool erase_flag;

        if (if_freq_mode_or) {
            erase_flag = true;
            for (auto &c : cohorts) {
                if (c.sumstat(iter->second, 3) >= freq_threshold && c.sumstat(iter->second, 3) <= 1-freq_threshold) {
                    erase_flag = false; 
                    break;
                }
            }
        } else {
            erase_flag = false;
            for (auto &c : cohorts) {
                if (c.sumstat(iter->second, 3) < freq_threshold || c.sumstat(iter->second, 3) > 1-freq_threshold) {
                    erase_flag = true; 
                    break;
                }
            }
        }
        
        if (erase_flag)
            iter = commonSNP_index_map.erase(iter);
        else {
            for (auto &c : cohorts)
                c.sumstat_screened.row(temp_index) = c.sumstat.row(iter->second);
            
            final_commonSNP.push_back(iter->first);
            final_commonSNP_index.push_back(iter->second);
            temp_index++;
            iter++;
        }
    }
    
    if (temp_index == 0)
        LOGGER.e(0, "Input data has no common SNPs");

    for (int n=0; n<cohort_num; n++) {
        auto &c = cohorts[n];
        c.sumstat_screened.conservativeResize(temp_index, 7);
        c.calc_Vp();
        c.sumstat = c.sumstat_screened;
        LOGGER << "Vp Cohort " << n+1 << ": " << c.Vp << endl;
    }
      
    LOGGER << endl;
    LOGGER.i(0, "common SNPs at last", to_string(temp_index));
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "--------------------------------" << endl << endl;
}
