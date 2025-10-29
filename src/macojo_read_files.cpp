#include "macojo.h"


void MACOJO::read_SNP_only(string filename, vector<string> &SNP_list, bool if_sumstat, bool if_bim) 
{
    ifstream sFile(filename.c_str());
    if (!sFile) LOGGER.e(0, "cannot open file [" + filename + "] to read");

    string str_buf, SNP_buf, temp;

    if (if_sumstat) getline(sFile, str_buf); // sumstat files have header line

    while (getline(sFile, str_buf)) {
        if (str_buf.empty()) continue; // skip empty lines

        istringstream iss(str_buf);
        if (if_bim) 
            iss >> temp >> SNP_buf; // bim file: chr, SNP, cm, pos, A1, A2
        else    
            iss >> SNP_buf; // assume the first column is SNP

        SNP_list.push_back(SNP_buf);
    }

    sFile.close();

    sort(SNP_list.begin(), SNP_list.end());
    if (adjacent_find(SNP_list.begin(), SNP_list.end()) != SNP_list.end())
        LOGGER.e(0, "Duplicate SNP in [" + filename + "], please check", *adjacent_find(SNP_list.begin(), SNP_list.end()));
    
    LOGGER << SNP_list.size() << " SNPs in file [" + filename + "]" << endl;
}


void Cohort::skim_fam(string PLINKfile) 
{   
    // read .fam to get individual number 
    string famFile = PLINKfile+".fam";
    ifstream Fam(famFile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open PLINK FAM file [" + famFile + "] to read");

    string FID, IID, str_buf;
    vector<string> indi_ids;

    while (Fam >> FID >> IID >> str_buf >> str_buf >> str_buf >> str_buf) {
        indi_ids.push_back(FID+':'+IID);
    }

    Fam.close();

    sort(indi_ids.begin(), indi_ids.end());
    if (adjacent_find(indi_ids.begin(), indi_ids.end()) != indi_ids.end())
        LOGGER.e(0, "Duplicate individual ID in FAM file [" + famFile + "], please check");

    indi_num = indi_ids.size();
    LOGGER << indi_num << " individuals in FAM file [" + famFile + "]" << endl;
}


void MACOJO::set_reference_from_bim()
{
    for (size_t n = 0; n < cohorts.size(); n++) {
        string bimFile = bfile_list[n]+".bim";
        ifstream Bim(bimFile.c_str());
        if (!Bim) LOGGER.e(0, "cannot open PLINK BIM file [" + bimFile + "] to read");
        LOGGER << "Reading PLINK BIM file from [" + bimFile + "] ..." << endl;

        int chr_buf, ibuf, ref_index;
        string SNP_buf, A1_buf = "0", A2_buf = "0", str_buf;
        if (n > 0) cohorts[n].swap_array.assign(shared.total_SNP_num, false);
        
        while (Bim >> chr_buf >> SNP_buf >> str_buf >> ibuf >> A1_buf >> A2_buf) {
            auto iter = shared.goodSNP_index_map.find(SNP_buf);
            if (iter == shared.goodSNP_index_map.end()) {
                cohorts[n].bim_index_list.push_back(-1);
                continue;
            }

            ref_index = iter->second;
            to_upper(A1_buf);
            to_upper(A2_buf);

            if (n == 0) {
                // first cohort, set as reference
                shared.chr_ref[ref_index] = chr_buf;
                shared.A1_ref[ref_index] = A1_buf;
                shared.A2_ref[ref_index] = A2_buf;
                shared.SNP_pos_ref[ref_index] = ibuf;
            } else {
                // other cohorts, check consistency
                if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
                else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {cohorts[n].swap_array[ref_index] = true;}
                else {
                    LOGGER.w(1, "removed, A1 and A2 different between two BIM files, please check", SNP_buf);
                    shared.goodSNP_index_map.erase(iter);
                    continue;
                }

                if (ibuf != shared.SNP_pos_ref[ref_index]) {
                    LOGGER.w(1, "removed, SNP position different between two BIM files, please check", SNP_buf);
                    shared.goodSNP_index_map.erase(iter);
                    continue;
                }      
            }

            cohorts[n].bim_index_list.push_back(ref_index);
        }

        Bim.close();
    }
}


void Cohort::read_sumstat(string sumstatfile) 
{  
    ifstream Meta(sumstatfile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open sumstat file [" + sumstatfile + "] to read");
    LOGGER << "Reading GWAS summary-level statistics from [" + sumstatfile + "] ..." << endl;

    string SNP_buf, A1_buf, A2_buf, freq_buf, b_buf, se_buf, p_buf, N_buf, str_buf;
    double freq, b, se, p, N, h;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line
    if (split_string(str_buf, vs_buf) != 8) LOGGER.e(0, "format error in sumstat file [" + sumstatfile + "]");
    
    map<string, int>::iterator iter;
    int ref_index;
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
            iter = shared.goodSNP_index_map.find(SNP_buf);
            if (iter != shared.goodSNP_index_map.end()) {
                LOGGER.w(1, "removed, invalid value in sumstat file [" + sumstatfile + "]", SNP_buf);
                shared.goodSNP_index_map.erase(iter);
            }
            continue;
        } 

        h = 2 * freq * (1 - freq);
        Vp_gcta_list.push_back(h * N * (se * se + b * b / (N - 1.0)));

        iter = shared.goodSNP_index_map.find(SNP_buf);
        if (iter == shared.goodSNP_index_map.end()) continue;

        ref_index = iter->second;

        if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
        else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {freq = 1 - freq; b = -b;}
        else {
            LOGGER.w(1, "removed, A1 and A2 different from ref BIM file, please check", SNP_buf);
            shared.goodSNP_index_map.erase(iter);
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

    Vp = median(Vp_gcta_list);
    sumstat.col(5) = sumstat.col(3) * (1 - sumstat.col(3)) * 2;
    sumstat.col(4) = (Vp / sumstat.col(5) - square(sumstat.col(0))) / sumstat.col(1) + 1;
    sumstat.col(6) = sumstat.col(4) * sumstat.col(5);
    LOGGER << "Estimated Vp: " << Vp << endl;
}


void Cohort::skim_frq(string PLINKfile) 
{   
    string frqFile = PLINKfile+".frq";
    ifstream Frq(frqFile.c_str());
    if (!Frq) {
        LOGGER.i(0, "cannot open PLINK FRQ file [" + frqFile + "] to read, skip allele frequency check");
        return;
    }

    LOGGER << "Reading PLINK FRQ file from [" + frqFile + "] ..." << endl;
    string SNP_buf, A1_buf, A2_buf, freq_buf, str_buf;
    double freq;

    vector<string> vs_buf;
    getline(Frq, str_buf); // skip header line
    if (split_string(str_buf, vs_buf) != 6) LOGGER.e(0, "format error in frq file [" + frqFile + "]");

    map<string, int>::iterator iter;
    int ref_index;

    while (Frq >> str_buf >> SNP_buf >> A1_buf >> A2_buf >> freq_buf >> str_buf) {
        to_upper(A1_buf);
        to_upper(A2_buf);
        freq = atof(freq_buf.c_str());

        iter = shared.goodSNP_index_map.find(SNP_buf);
        if (iter == shared.goodSNP_index_map.end()) continue;

        ref_index = iter->second;

        if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
        else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {freq = 1 - freq;}
        else {
            LOGGER.w(1, "removed, A1 and A2 different from ref BIM file, please check", SNP_buf);
            shared.goodSNP_index_map.erase(iter);
            continue;
        } 

        if (abs(sumstat(ref_index, 3) - freq) > params.diff_freq) {
            // LOGGER.w(1, "removed, allele frequency too different between sumstat and bedfile", SNP_buf);
            shared.goodSNP_index_map.erase(SNP_buf);
            continue;
        }
    }

    Frq.close();
}

void Cohort::read_PLINK_bed(string PLINKfile, bool is_ref_cohort) 
{   
    string bedFile = PLINKfile+".bed";
    ifstream Bed(bedFile.c_str(), ios::binary | ios::ate);
    if (!Bed) LOGGER.e(0, "cannot open PLINK BED file [" + bedFile + "] to read");    
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    
    const int bytes_per_snp = (indi_num + 3) / 4;

    if (Bed.tellg() != streamoff(uint64_t(bytes_per_snp) * bim_index_list.size() + 3)) 
        LOGGER.e(0, "PLINK BED file [" + bedFile + "] size does not match .bim and .fam files, please check");

    Bed.seekg(0, ios::beg);
    char ch[3];
    for (int i=0; i<3; i++) {Bed.read(&ch[i], 1);} 
    if (!Bed || ch[0] != 0x6C || ch[1] != 0x1B || ch[2] != 0x01)
        LOGGER.e(0, "PLINK BED file [" + bedFile + "] not in SNP-major mode, please check");

    vector<char> buffer(bytes_per_snp);

    for (int& ref_index : bim_index_list) {
        Bed.read(buffer.data(), bytes_per_snp);
        if (ref_index == -1) continue; // SNP not in common SNP list

        bool swap = is_ref_cohort ? false : swap_array[ref_index];
        string SNP_buf = shared.SNP_ref[ref_index];

        if (!genotype.decode_single_genotype(buffer, ref_index, swap)) {
            LOGGER.w(1, "removed, all values are NA or identical in bedfile", SNP_buf);
            shared.goodSNP_index_map.erase(SNP_buf);
            ref_index = -1;
            continue;
        }
        
        if (fabs(sumstat(ref_index, 3) - genotype.X_avg[ref_index]/2) > params.diff_freq) {
            // LOGGER.w(1, "removed, allele frequency too different between sumstat and bedfile", SNP_buf);
            shared.goodSNP_index_map.erase(SNP_buf);
            ref_index = -1;
            continue;
        }

        if (genotype.missingness[ref_index] > params.missingness) {
            // LOGGER.w(1, "removed, missingness too high in bedfile", SNP_buf);
            shared.goodSNP_index_map.erase(SNP_buf);
            ref_index = -1;
            continue;
        }
    }

    Bed.close();
    LOGGER << "Finished reading PLINK BED file" << endl;
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

    while (Ld >> SNP1_buf >> SNP2_buf >> r_buf) {

        if ((iter = shared.goodSNP_index_map.find(SNP1_buf)) == shared.goodSNP_index_map.end()) continue;
        SNP1_index = iter->second;

        if ((iter = shared.goodSNP_index_map.find(SNP2_buf)) == shared.goodSNP_index_map.end()) continue;
        SNP2_index = iter->second;

        LD_matrix(SNP1_index, SNP2_index) = atof(r_buf.c_str()); 
    }

    Ld.close();

    if (!is_ref_cohort) {
        // adjust LD_packed according to allele coding
        for (int i = 0; i < shared.total_SNP_num; i++) {
            for (int j = i; j < shared.total_SNP_num; j++) {
                if (swap_array[i] ^ swap_array[j]) {
                    LD_matrix(i,j) = -LD_matrix(i,j);
                    cout << "swap " << i << " " << j << endl;
                }
            }
        }
    }

    LOGGER << "Finished reading PLINK LD file" << endl;
}


void MACOJO::read_cojo_PLINK_files() 
{   
    int cohort_num = cojo_file_list.size();
    vector<string> commonSNP_all_cohorts;

    // if user provides extract SNP file, read file and initialize common SNPs
    if (!extract_file.empty()) {
        read_SNP_only(extract_file, commonSNP_all_cohorts);
        if (commonSNP_all_cohorts.size() == 0)
            LOGGER.e(0, "Input extract SNP file has no SNPs", extract_file);

        LOGGER.i(0, "SNPs initialized by the user\n", to_string(commonSNP_all_cohorts.size()));
    }

    // Step 1: get common SNPs across all cohorts
    for (int n = 0; n < cohort_num; n++) {
        auto start = chrono::steady_clock::now();

        vector<string> SNP_sumstat, SNP_PLINK, temp;

        string sumstat_file = cojo_file_list[n], PLINK_file = bfile_list[n];
        read_SNP_only(PLINK_file+".bim", SNP_PLINK, false, true);

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
        
        auto end = chrono::steady_clock::now();
        LOGGER.i(0, "common SNPs after including Cohort "+to_string(n+1), to_string(commonSNP_all_cohorts.size()));
        LOGGER << "Time taken: " << chrono::duration<double>(end-start).count() << " seconds" << endl << endl;
    }

    if (commonSNP_all_cohorts.size() == 0)
        LOGGER.e(0, "Input data has no common SNPs across all cohorts");
    
    // Step 2: set reference based on BIM file of cohort 1
    auto start = chrono::steady_clock::now();

    int temp_index = 0;
    for (const auto& snp : commonSNP_all_cohorts) {
        shared.goodSNP_index_map.insert(make_pair(snp, temp_index));
        temp_index++;
    }
    
    shared.SNP_ref = commonSNP_all_cohorts;
    shared.total_SNP_num = temp_index;
    shared.chr_ref.resize(temp_index);
    shared.A1_ref.resize(temp_index);
    shared.A2_ref.resize(temp_index);
    shared.SNP_pos_ref.resize(temp_index);

    set_reference_from_bim();
    
    auto end = chrono::steady_clock::now();
    LOGGER << "Time taken: " << chrono::duration<double>(end-start).count() << " seconds" << endl << endl;

    start = chrono::steady_clock::now();
    // Step 3: read sumstat files for all cohorts
    for (int n = 0; n < cohort_num; n++) {
        auto& c = cohorts[n];
        // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D
        c.sumstat.resize(shared.total_SNP_num, 7);
        c.read_sumstat(cojo_file_list[n]);
    }

    // exclude rare SNPs based on freq_threshold
    for (auto iter = shared.goodSNP_index_map.begin(); iter != shared.goodSNP_index_map.end(); ) {
        int ref_index = iter->second;
        bool erase_flag = !params.if_freq_mode_and;

        if (params.if_freq_mode_and) {
            for (auto &c : cohorts) {
                if (c.sumstat(ref_index, 3) < params.maf || c.sumstat(ref_index, 3) > 1-params.maf) {
                    erase_flag = true; 
                    break;
                }
            }
        } else {
            for (auto &c : cohorts) {
                if (c.sumstat(ref_index, 3) >= params.maf && c.sumstat(ref_index, 3) <= 1-params.maf) {
                    erase_flag = false; 
                    break;
                }
            }
        } 
        
        if (erase_flag) {
            LOGGER.w(1, "removed, rare SNP based on frequency threshold", shared.SNP_ref[ref_index]);
            iter = shared.goodSNP_index_map.erase(iter);
        }
        else 
            iter++;
    }

    LOGGER.i(0, "common SNPs after removing rare SNPs", to_string(shared.goodSNP_index_map.size()));

    end = chrono::steady_clock::now();
    LOGGER << "Time taken: " << chrono::duration<double>(end-start).count() << " seconds" << endl << endl;

    if (shared.goodSNP_index_map.size() == 0)
        LOGGER.e(0, "Input data has no common SNPs");
    
    // Step 4: read bim and bed files for all cohorts
    for (int n = 0; n < cohort_num; n++) {
        // bim file of cohort 1 is reference, do not need to check allele
        auto start = chrono::steady_clock::now();
        auto& c = cohorts[n];
    
        if (params.if_LD_mode) {
            c.skim_frq(bfile_list[n]);
            c.LD_matrix.resize(shared.total_SNP_num);
            c.read_PLINK_LD(bfile_list[n], n==0);
        }
        else {
            c.skim_fam(bfile_list[n]);
            c.genotype.resize(shared.total_SNP_num, c.indi_num);
            c.read_PLINK_bed(bfile_list[n], n==0);
        }

        auto end = chrono::steady_clock::now();
        LOGGER.i(0, "common SNPs remaining", to_string(shared.goodSNP_index_map.size()));
        LOGGER << "Time taken: " << chrono::duration<double>(end-start).count() << " seconds" << endl << endl;

        if (shared.goodSNP_index_map.size() == 0)
            LOGGER.e(0, "Input data has no common SNPs");
    }

    // Step 5: finalize SNP index lists
    for (auto iter = shared.goodSNP_index_map.begin(); iter != shared.goodSNP_index_map.end(); iter++) 
        good_SNP.push_back(iter->second);

    vector<int> count_array(shared.total_SNP_num);
    for (int i = 0; i < shared.total_SNP_num; i++)
        count_array[i] = i;
    
    set_difference(count_array.begin(), count_array.end(), 
        good_SNP.begin(), good_SNP.end(), back_inserter(bad_SNP));

    screened_SNP = good_SNP;

    LOGGER.i(0, "common SNPs at last for iterative selection", to_string(shared.goodSNP_index_map.size()));
    LOGGER << "--------------------------------" << endl << endl;
}
