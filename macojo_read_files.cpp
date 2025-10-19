#include "macojo.h"


namespace {
    struct Entry {
        bool A1, A2;
        uint8_t geno;
        bool valid;
    };

    Entry LUT_normal[256][4];
    Entry LUT_swapped[256][4];

    void build_LUT() 
    {
        for (int b = 0; b < 256; b++) {
            for (int k = 0; k < 4; k++) {
                int code = (b >> (2*k)) & 0x3;
                Entry e;
                switch(code) {
                    case 0: e = {1,1,2,true}; break;  // 00 hom major
                    case 2: e = {0,1,1,true}; break;  // 10 hetero
                    case 3: e = {0,0,0,true}; break;  // 11 hom minor
                    case 1: e = {1,0,0,false}; break; // 01 missing
                }
                LUT_normal[b][k] = e;

                // swap: 00<->11
                Entry e2 = e;
                if (code == 0) e2 = {0,0,0,true};
                else if (code == 3) e2 = {1,1,2,true};
                LUT_swapped[b][k] = e2;
            }
        }
    }
};


void MACOJO::read_SNP_only(string filename, vector<string> &SNP_list, bool if_sumstat, bool if_bim) 
{
    ifstream sFile(filename.c_str());
    if (!sFile) LOGGER.e(0, "cannot open file [" + filename + "] to read");

    string str_buf, SNP_buf, temp;
    int SNP_num = 0;

    if (if_sumstat) getline(sFile, str_buf); // sumstat files have header line

    while (getline(sFile, str_buf)) {
        if (str_buf.empty()) continue; // skip empty lines

        istringstream iss(str_buf);
        if (if_bim) 
            iss >> temp >> SNP_buf; // bim file: chr, SNP, cm, pos, A1, A2
        else    
            iss >> SNP_buf; // assume the first column is SNP

        SNP_list.push_back(SNP_buf);
        SNP_num++;
    }

    sFile.close();

    sort(SNP_list.begin(), SNP_list.end());
    if (adjacent_find(SNP_list.begin(), SNP_list.end()) != SNP_list.end())
        LOGGER.e(0, "Duplicate SNP in [" + filename + "], please check", *adjacent_find(SNP_list.begin(), SNP_list.end()));
    
    LOGGER << SNP_num << " SNPs in file [" + filename + "]" << endl;
}


void Cohort::skim_fam(string famFile) 
{   
    // read .fam to get individual number 
    ifstream Fam(famFile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open FAM file [" + famFile + "] to read");

    string FID, IID, str_buf;
    vector<string> indi_ids;
    indi_num = 0;

    while (Fam >> FID >> IID >> str_buf >> str_buf >> str_buf >> str_buf) {
        indi_ids.push_back(FID+':'+IID);
        indi_num++;
    }

    Fam.close();

    sort(indi_ids.begin(), indi_ids.end());
    if (adjacent_find(indi_ids.begin(), indi_ids.end()) != indi_ids.end())
        LOGGER.e(0, "Duplicate individual ID in FAM file [" + famFile + "], please check");

    LOGGER << indi_num << " individuals in FAM file [" + famFile + "]" << endl;
}


void MACOJO::set_reference_from_bim(string PLINKfile)
{
    string bimFile = PLINKfile+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "] to read");

    // Read bim file
    int chr_buf, ibuf, ref_index;
    string SNP_buf, A1_buf = "0", A2_buf = "0", str_buf;
    map<string, int>::iterator iter;
    
    while (Bim >> chr_buf >> SNP_buf >> str_buf >> ibuf >> A1_buf >> A2_buf) {
        iter = shared.commonSNP_index_map.find(SNP_buf);
        if (iter == shared.commonSNP_index_map.end()) continue;

        to_upper(A1_buf);
        to_upper(A2_buf);

        ref_index = iter->second;
        shared.chr_ref[ref_index] = chr_buf;
        shared.A1_ref[ref_index] = A1_buf;
        shared.A2_ref[ref_index] = A2_buf;
        shared.SNP_pos_ref[ref_index] = ibuf;
    }

    Bim.close();
}


void Cohort::read_sumstat(string sumstatfile) 
{  
    ifstream Meta(sumstatfile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open the file [" + sumstatfile + "] to read");

    string SNP_buf, A1_buf, A2_buf, freq_buf, b_buf, se_buf, p_buf, N_buf, str_buf;
    double freq, b, se, p, N, h;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line
    if (split_string(str_buf, vs_buf) < 7) LOGGER.e(0, "format error in sumstat file [" + sumstatfile + "]");
    
    map<string, int>::iterator iter;
    int ref_index;

    sumstat_screened.resize(shared.commonSNP_total_num, 9);
    sumstat_screened.col(8).setOnes(); // col 8: X_norm_square, initialized as 1
    vector<double> Vp_gcta_list;

    while (Meta >> SNP_buf >> A1_buf >> A2_buf >> freq_buf >> b_buf >> se_buf >> p_buf >> N_buf) {

        to_upper(A1_buf);
        to_upper(A2_buf);
        freq = atof(freq_buf.c_str());
        b = atof(b_buf.c_str());
        se = atof(se_buf.c_str());
        p = atof(p_buf.c_str());
        N = atof(N_buf.c_str());
        
        if ((b_buf == "NA" || b_buf == "." || se_buf == "NA" || se_buf == "." || se_buf == "0" || \
            p_buf == "NA" || p_buf == "." || N_buf == "NA" || N_buf == "." || N < 10)) { 
            iter = shared.commonSNP_index_map.find(SNP_buf);
            if (iter != shared.commonSNP_index_map.end()) {
                LOGGER.w(1, "removed, invalid value in sumstat file [" + sumstatfile + "]", SNP_buf);
                shared.commonSNP_index_map.erase(iter);
            }
            continue;
        } 

        h = 2 * freq * (1 - freq);
        Vp_gcta_list.push_back(h * N * (se * se + b * b / (N - 1.0)));

        iter = shared.commonSNP_index_map.find(SNP_buf);
        if (iter == shared.commonSNP_index_map.end()) continue;

        ref_index = iter->second;

        if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
        else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {freq = 1 - freq; b = -b;}
        else {
            LOGGER.w(1, "removed, A1 and A2 different from ref BIM file, please check", SNP_buf);
            shared.commonSNP_index_map.erase(iter);
            continue;
        } 

        // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D, 7:X_avg, 8:X_norm_square
        sumstat_screened(ref_index, 0) = b;
        sumstat_screened(ref_index, 1) = se * se;
        sumstat_screened(ref_index, 2) = p;
        sumstat_screened(ref_index, 3) = freq;
        sumstat_screened(ref_index, 4) = N;
    }

    Meta.close();
    Vp = median(Vp_gcta_list);
}


void Cohort::read_PLINK(string PLINKfile, bool is_ref_cohort) 
{   
    string bimFile = PLINKfile+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "] to read");
    LOGGER << "Reading PLINK BIM file from [" + bimFile + "]..." << endl;

    int ibuf = 0;
    string SNP_buf, A1_buf = "0", A2_buf = "0", str_buf;

    // check bed file header
    string bedFile = PLINKfile+".bed";
    fstream Bed(bedFile.c_str(), ios::in | ios::binary);
    if (!Bed) LOGGER.e(0, "cannot open BED file [" + bedFile + "] to read");

    char ch[3];
    for (int i=0; i<3; i++) {Bed.read(&ch[i], 1);} 
    if (!Bed || ch[0] != 0x6C || ch[1] != 0x1B || ch[2] != 0x01)
        LOGGER.e(0, "PLINK BED file [" + bedFile + "] not in SNP-major mode, please check");
        
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;

    const uint64_t words_per_snp = (indi_num + 63) / 64;
    const int bytes_per_snp = (indi_num + 3) / 4;
    const int full_bytes = indi_num / 4;
    const int remainder  = indi_num % 4;

    X_A1.resize(shared.commonSNP_total_num * words_per_snp);
    X_A2.resize(shared.commonSNP_total_num * words_per_snp);
    vector<unsigned char> buffer(bytes_per_snp);

    int bad_freq_count = 0;

    while (Bim >> str_buf >> SNP_buf >> str_buf >> ibuf >> A1_buf >> A2_buf) {
  
        Bed.read(reinterpret_cast<char*>(buffer.data()), bytes_per_snp);
        if (!Bed) break;

        auto iter = shared.commonSNP_index_map.find(SNP_buf);
        if (iter == shared.commonSNP_index_map.end()) continue;

        int ref_index = iter->second;
        bool swap = false;
        
        // check allele and position in BIM file with reference cohort
        if (!is_ref_cohort) {
            to_upper(A1_buf);
            to_upper(A2_buf);

            if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
            else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {swap = true;}
            else {
                LOGGER.w(1, "removed, A1 and A2 different between two BIM files, please check", SNP_buf);
                shared.commonSNP_index_map.erase(iter);
                continue;
            }

            if (ibuf != shared.SNP_pos_ref[ref_index]) {
                LOGGER.w(1, "removed, SNP position different between two BIM files, please check", SNP_buf);
                shared.commonSNP_index_map.erase(iter);
                continue;
            }      
        }
        
        // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
        const auto& LUT = swap ? LUT_swapped : LUT_normal;
        const uint64_t base_word_index = uint64_t(ref_index) * words_per_snp;

        double SNP_sum = 0, SNP_square_sum = 0, not_NA_indi_num = 0;
        uint64_t wordA1 = 0ULL, wordA2 = 0ULL, word_index = 0;

        int bit_index  = 0;
        for (int j = 0; j < full_bytes; j++) {
            unsigned char b = buffer[j];
            const auto &lut = LUT[b];

            // unrolled 4 genotypes
            for (int k = 0; k < 4; k++) {
                const auto &e = lut[k];
                const uint64_t mask = 1ULL << (bit_index & 63);
                wordA1 |= uint64_t(e.A1) * mask;
                wordA2 |= uint64_t(e.A2) * mask;

                const double g = e.geno;
                const double v = double(e.valid);
                SNP_sum += g * v;
                SNP_square_sum += g * g * v;
                not_NA_indi_num += v;

                bit_index++;
                if ((bit_index & 63ULL) == 0ULL) {
                    X_A1[base_word_index + word_index] = wordA1;
                    X_A2[base_word_index + word_index] = wordA2;
                    wordA1 = wordA2 = 0ULL;
                    word_index++;
                }
            }        
        }

        // Tail: remaining 1–3 individuals
        if (remainder > 0) {
            unsigned char b = buffer[full_bytes];
            const auto &lut = LUT[b];
            for (int k = 0; k < remainder; k++) {
                const auto &e = lut[k];
                const uint64_t mask = 1ULL << (bit_index & 63);
                wordA1 |= uint64_t(e.A1) * mask;
                wordA2 |= uint64_t(e.A2) * mask;

                const double g = e.geno;
                const double v = double(e.valid);
                SNP_sum += g * v;
                SNP_square_sum += g * g * v;
                not_NA_indi_num += v;
                bit_index++;
            }
        }

        // handle leftover partial word
        if ((bit_index & 63ULL) != 0ULL) {
            X_A1[base_word_index + word_index] = wordA1;
            X_A2[base_word_index + word_index] = wordA2;
        }

        if (not_NA_indi_num == 0) {
            LOGGER.w(1, "removed, all values are NA in bedfile", SNP_buf);
            shared.commonSNP_index_map.erase(iter);
            continue;
        }

        double SNP_avg = SNP_sum/not_NA_indi_num;
        if (fabs(sumstat_screened(ref_index, 3) - SNP_avg/2) > params.freq_diff_threshold) {
            // LOGGER.w(1, "removed, allele frequency too different between sumstat and bedfile", SNP_buf);
            bad_freq_count++;
            shared.commonSNP_index_map.erase(iter);
            continue;
        }
        
        double SSPD = not_NA_indi_num*SNP_square_sum - SNP_sum*SNP_sum;
        if (SSPD < 0.5) {
            LOGGER.w(1, "removed, identical genotypes for all individuals in bedfile", SNP_buf);
            shared.commonSNP_index_map.erase(iter);
            continue;
        }
        
        if (!params.if_keep_NA) {
            sumstat_screened(ref_index, 7) = SNP_avg;
            sumstat_screened(ref_index, 8) = SSPD/not_NA_indi_num; // this is actually variance*(n-1)
        }
    }

    Bim.close();
    Bed.close();

    LOGGER << "Finished reading PLINK file" << endl;
    if (bad_freq_count > 0)
        LOGGER.w(1, "SNPs removed due to allele frequency too different between sumstat and bedfile", to_string(bad_freq_count));
}


void Cohort::read_PLINK_LD(string PLINKfile, bool is_ref_cohort) 
{   
    string ldFile = PLINKfile+".ld";
    ifstream Ld(ldFile.c_str());
    if (!Ld) LOGGER.e(0, "cannot open PLINK LD file [" + ldFile + "] to read");
    LOGGER << "Reading PLINK LD file from [" + ldFile + "] ..." << endl;

    // Read ld file
    string SNP1_buf, SNP2_buf, r_buf;
    int SNP1_index, SNP2_index;
    map<string, int>::iterator iter;

    LD_matrix.resize(shared.commonSNP_total_num);

    while (Ld >> SNP1_buf >> SNP2_buf >> r_buf) {

        if ((iter = shared.commonSNP_index_map.find(SNP1_buf)) == shared.commonSNP_index_map.end()) continue;
        SNP1_index = iter->second;

        if ((iter = shared.commonSNP_index_map.find(SNP2_buf)) == shared.commonSNP_index_map.end()) continue;
        SNP2_index = iter->second;

        LD_matrix(SNP1_index, SNP2_index) = atof(r_buf.c_str()); 
    }

    Ld.close();

    if (!is_ref_cohort) {
        string bimFile = PLINKfile+".bim";
        ifstream Bim(bimFile.c_str());
        if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "] to read");
        LOGGER << "Adjusting LD values according to PLINK BIM file from [" + bimFile + "]..." << endl;
        
        // Read bim file
        int ibuf = 0;
        string SNP_buf, A1_buf = "0", A2_buf = "0", str_buf;
    
        vector<int> swap_array(shared.commonSNP_total_num, 0);

        while (Bim >> str_buf >> SNP_buf >> str_buf >> ibuf >> A1_buf >> A2_buf) {

            if ((iter = shared.commonSNP_index_map.find(SNP_buf)) == shared.commonSNP_index_map.end()) continue;
            
            int ref_index = iter->second;

            // check allele and position in BIM file with reference cohort
            to_upper(A1_buf);
            to_upper(A2_buf);

            if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
            else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {swap_array[ref_index] = 1;}
            else {
                LOGGER.w(1, "removed, A1 and A2 different between two BIM files, please check", SNP_buf);
                shared.commonSNP_index_map.erase(iter);
                continue;
            }

            if (ibuf != shared.SNP_pos_ref[ref_index]) {
                LOGGER.w(1, "removed, SNP position different between two BIM files, please check", SNP_buf);
                shared.commonSNP_index_map.erase(iter);
                continue;
            }
        }
            
        Bim.close();
    
        // adjust LD_packed according to allele coding
        for (int i = 0; i < shared.commonSNP_total_num; i++) {
            for (int j = i; j < shared.commonSNP_total_num; j++) {
                if (swap_array[i] ^ swap_array[j]) {
                    LD_matrix(i,j) = -LD_matrix(i,j);
                    cout << "swap " << i << " " << j << endl;
                }
            }
        }
    }

    LOGGER << "Finished reading PLINK LD file" << endl;
}


void Cohort::calc_adjusted_N() 
{   
    ArrayXXd &s = sumstat_screened;

    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D, 7:X_avg, 8:X_norm_square 
    s.col(5) = s.col(3) * (1-s.col(3)) * 2;

    if (!params.if_gcta_COJO) {
        ArrayXd Vp_gcta_list = s.col(5) * s.col(4) * (s.col(1) + square(s.col(0))/(s.col(4)-1));
        Vp = median(Vp_gcta_list);
    }

    s.col(4) = (Vp / s.col(5) - square(s.col(0))) / s.col(1) + 1;
    s.col(6) = s.col(4) * s.col(5);

    if (!params.if_gcta_COJO) {
        ArrayXd Vp_gcta_list = s.col(6) * (s.col(1) + square(s.col(0))/(s.col(4)-1));
        Vp = median(Vp_gcta_list);
    }
}


void MACOJO::read_cojo_PLINK_files(char** filenames, int cohort_num, string extract_file) 
{   
    vector<string> commonSNP_all_cohorts;

    // if user provides extract SNP file, read file and initialize common SNPs
    if (!extract_file.empty()) {
        read_SNP_only(extract_file, commonSNP_all_cohorts);
        if (commonSNP_all_cohorts.size() == 0)
            LOGGER.e(0, "Input extract SNP file has no SNPs", extract_file);

        LOGGER.i(0, "SNPs initialized by the user\n", to_string(commonSNP_all_cohorts.size()));
    }

    // Step 1: get common SNPs across all cohorts
    clock_t tStart;

    for (int n=0; n<cohort_num; n++) {

        tStart = clock();

        auto &c = cohorts[n];
        vector<string> SNP_sumstat, SNP_PLINK, temp;

        string sumstat_file = filenames[n*2], PLINK_file = filenames[n*2+1];
        read_SNP_only(PLINK_file+".bim", SNP_PLINK, false, true);
        c.skim_fam(PLINK_file+".fam");

        if (n == 0 && commonSNP_all_cohorts.size() == 0)
            commonSNP_all_cohorts = SNP_PLINK;
        else {
            set_intersection(commonSNP_all_cohorts.begin(), commonSNP_all_cohorts.end(), 
                SNP_PLINK.begin(), SNP_PLINK.end(), back_inserter(temp));
            commonSNP_all_cohorts.swap(temp);
        }

        vector<string>().swap(temp); // clear temp
        read_SNP_only(sumstat_file, SNP_sumstat, true, false);
        set_intersection(commonSNP_all_cohorts.begin(), commonSNP_all_cohorts.end(), 
            SNP_sumstat.begin(), SNP_sumstat.end(), back_inserter(temp));
        commonSNP_all_cohorts.swap(temp);
        
        LOGGER.i(0, "common SNPs after including Cohort "+to_string(n+1), to_string(commonSNP_all_cohorts.size()));
        LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl << endl;
    } 

    if (commonSNP_all_cohorts.size() == 0)
        LOGGER.e(0, "Input data has no common SNPs across all cohorts");
    
        // Step 2: set reference based on BIM file of cohort 1
    tStart = clock();

    int temp_index = 0;

    for (const auto &snp : commonSNP_all_cohorts) {
        shared.commonSNP_index_map.insert(make_pair(snp, temp_index));
        temp_index++;
    }

    shared.chr_ref.resize(temp_index);
    shared.A1_ref.resize(temp_index);
    shared.A2_ref.resize(temp_index);
    shared.SNP_pos_ref.resize(temp_index);
    shared.commonSNP_total_num = temp_index;

    LOGGER << "Setting reference based on BIM file of Cohort 1 ..." << endl;
    set_reference_from_bim(filenames[1]);

    // Step 3: read sumstat files for all cohorts
    for (int n=0; n<cohort_num; n++) {
        string sumstat_file = filenames[n*2];
        LOGGER << "Reading Cohort " << n+1 << " GWAS summary-level statistics from [" + sumstat_file + "] ..." << endl;
        cohorts[n].read_sumstat(sumstat_file);
    }

    // exclude rare SNPs based on freq_threshold
    for (auto iter = shared.commonSNP_index_map.begin(); iter != shared.commonSNP_index_map.end(); ) {
        bool erase_flag = !params.if_freq_mode_and;

        if (params.if_freq_mode_and) {
            for (auto &c : cohorts) {
                if (c.sumstat_screened(iter->second, 3) < params.freq_threshold || 
                    c.sumstat_screened(iter->second, 3) > 1-params.freq_threshold) {
                    erase_flag = true; 
                    break;
                }
            }
        } else {
            for (auto &c : cohorts) {
                if (c.sumstat_screened(iter->second, 3) >= params.freq_threshold && 
                    c.sumstat_screened(iter->second, 3) <= 1-params.freq_threshold) {
                    erase_flag = false; 
                    break;
                }
            }
        } 
        
        if (erase_flag)
            iter = shared.commonSNP_index_map.erase(iter);
        else 
            iter++;
    }

    LOGGER.i(0, "common SNPs after removing rare SNPs", to_string(shared.commonSNP_index_map.size()));
    LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl << endl;

    if (shared.commonSNP_index_map.size() == 0)
        LOGGER.e(0, "Input data has no common SNPs");
    
    // Step 4: read bim and bed files for all cohorts
    build_LUT();

    for (int n=0; n<cohort_num; n++) {
        // bim file of cohort 1 is reference, do not need to check allele
        tStart = clock();
    
        if (params.if_LD_mode)
            cohorts[n].read_PLINK_LD(filenames[n*2+1], n==0);
        else {
            cohorts[n].read_PLINK(filenames[n*2+1], n==0);
            cohorts[n].X_avg = cohorts[n].sumstat_screened.col(7);
            cohorts[n].X_square = cohorts[n].sumstat_screened.col(8);
        }
           
        LOGGER << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl << endl;
    }

    shared.commonSNP_total_num = shared.commonSNP_index_map.size();    
    if (shared.commonSNP_index_map.size() == 0)
        LOGGER.e(0, "Input data has no qualified SNPs");

    LOGGER.i(0, "common SNPs at last for iterative selection", to_string(shared.commonSNP_total_num));

    // Step 5: clean sumstat for all cohorts and finalize common SNP list
    for (auto iter = shared.commonSNP_index_map.begin(); iter != shared.commonSNP_index_map.end(); iter++) {
        shared.final_commonSNP.push_back(iter->first);
        shared.final_commonSNP_index.push_back(iter->second);
    }
    
    ArrayXXd temp;
    for (int n=0; n<cohort_num; n++) {
        auto &c = cohorts[n];
        temp.setZero(shared.commonSNP_total_num, 9);

        #pragma omp parallel for
        for (int i = 0; i < shared.commonSNP_total_num; i++)
            temp.row(i) = c.sumstat_screened.row(shared.final_commonSNP_index[i]);

        c.sumstat_screened = temp; // 0,1,2,3,4,7,8 are filled, 5,6 are to be calculated
        c.calc_adjusted_N(); // 5,6 are filled, 4 is adjusted
        LOGGER << "Vp Cohort " << n+1 << ": " << c.Vp << endl;
    }

    LOGGER << "--------------------------------" << endl << endl;
    /*
    // check if X and LD give same results
    for (auto &c : cohorts) {
        for (int i = 0; i < commonSNP_index_map.size(); i++) {
            for (int j = i; j < commonSNP_index_map.size(); j++) {

                VectorXd X1_vec, X2_vec;

                c.get_vector_from_bed_matrix(final_commonSNP_index[i], X1_vec);
                c.get_vector_from_bed_matrix(final_commonSNP_index[j], X2_vec);
                double r_X1_X2 = (X1_vec.transpose() * X2_vec).value() / (c.indi_num-1);

                double r_LD = c.LD_matrix(final_commonSNP_index[i], final_commonSNP_index[j]);

                if (fabs(r_X1_X2 - r_LD) > 1e-10) 
                    LOGGER.e(0, "LD and X give different results, please check");
            }
        }
    }
    */
}
