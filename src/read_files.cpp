#include "macojo.h"


void skim_fam(string filename, vector<string> &str_list) 
{   
    ifstream sFile(filename.c_str());
    if (!sFile) LOGGER.e("Cannot open file [" + filename + "] to read");

    string str_buf, FID, IID;

    while (getline(sFile, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        parse_string(pt, FID);
        parse_string(pt, IID);
        str_list.push_back(FID + ":" + IID);
    }

    sFile.close();

    sort(str_list.begin(), str_list.end());
    if (adjacent_find(str_list.begin(), str_list.end()) != str_list.end())
        LOGGER.e("Duplicate names in [" + filename + "], please check", *adjacent_find(str_list.begin(), str_list.end()));
}


// start from 1 for col_idx
void skim_file(string filename, vector<string> &str_list, int col_idx, bool has_header)
{
    ifstream sFile(filename.c_str());
    if (!sFile) LOGGER.e("Cannot open file [" + filename + "] to read");

    string str_buf;
    if (has_header) getline(sFile, str_buf); // skip header line

    while (getline(sFile, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        
        // iterate to column col_idx
        int col = 1;
        while (col < col_idx) {
            skip_token(pt);
            skip_delim(pt);
            if (!*pt) LOGGER.e("Requested column exceeds the total number of columns", filename);
            col++;
        }

        // p now at start of desired column
        const char* start = pt;
        skip_token(pt);
        str_list.emplace_back(start, pt - start);
    }

    sFile.close();

    sort(str_list.begin(), str_list.end());
    if (adjacent_find(str_list.begin(), str_list.end()) != str_list.end())
        LOGGER.e("Duplicate names in [" + filename + "], please check", *adjacent_find(str_list.begin(), str_list.end()));
}


void Cohort::read_sumstat() 
{  
    string cojo_file = params.cojo_file_list[cohort_index];
    ifstream Meta(cojo_file.c_str());
    if (!Meta) LOGGER.e("Cannot open sumstat file [" + cojo_file + "] to read");
    LOGGER.i("Reading GWAS summary-level statistics from [" + cojo_file + "] ...");

    string str_buf, SNP_buf, A1_buf, A2_buf;
    double freq, b, se, p, N;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line

    istringstream iss(str_buf);
    string token;
    while (iss >> token) vs_buf.push_back(token);
    if (vs_buf.size() != 8 || vs_buf[0] != "SNP" || vs_buf[1] != "A1" || vs_buf[2] != "A2" ||
        vs_buf[3] != "freq" || vs_buf[4] != "b" || vs_buf[5] != "se" || vs_buf[6] != "p" || vs_buf[7] != "N")
        LOGGER.e("Format error in sumstat file, please check");
    
    map<string, int>::iterator iter;
    int ref_index;
    vector<double> Vp_gcta_list;

    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D
    sumstat.resize(shared.total_SNP_num, 7);
    
    while (getline(Meta, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        parse_string(pt, SNP_buf);
        parse_string(pt, A1_buf, true);
        parse_string(pt, A2_buf, true);

        if (!parse_double(pt, freq) || !parse_double(pt, b) || !parse_double(pt, se) || !parse_double(pt, p) || !parse_double(pt, N)) { 
            iter = shared.goodSNP_index_map.find(SNP_buf);
            if (iter != shared.goodSNP_index_map.end()) {
                LOGGER.w("removed, invalid value in sumstat file", SNP_buf);
                shared.goodSNP_index_map.erase(iter);
            }
            continue;
        } 

        Vp_gcta_list.push_back(2 * freq * (1 - freq) * N * (se * se + b * b / (N - 1.0)));

        auto iter = shared.goodSNP_index_map.find(SNP_buf);
        if (iter == shared.goodSNP_index_map.end()) continue;

        ref_index = iter->second;

        if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
        else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {freq = 1 - freq; b = -b;}
        else {
            LOGGER.w("removed, A1 and A2 different from ref BIM file, please check", SNP_buf);
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


void Cohort::read_fam() 
{   
    string famFile = params.bfile_list[cohort_index]+".fam";
    LOGGER.i("Reading PLINK FAM file from [" + famFile + "] ...");

    vector<string> all_ids;
    skim_fam(famFile, all_ids);
    fam_indi_num = all_ids.size();

    if (!params.keep_file_list[cohort_index].empty()) {
        vector<string> keep_ids, temp_ids;
        skim_fam(params.keep_file_list[cohort_index], keep_ids);
        set_intersection(all_ids.begin(), all_ids.end(), keep_ids.begin(), keep_ids.end(), back_inserter(temp_ids));
        all_ids = temp_ids;
        LOGGER.i("individuals after keeping individuals", to_string(all_ids.size()));
    }

    if (!params.remove_file_list[cohort_index].empty()) {
        vector<string> remove_ids, temp_ids;
        skim_fam(params.remove_file_list[cohort_index], remove_ids);
        set_difference(all_ids.begin(), all_ids.end(), remove_ids.begin(), remove_ids.end(), back_inserter(temp_ids));
        all_ids = temp_ids;
        LOGGER.i("individuals after removing individuals", to_string(all_ids.size()));
    }

    if (all_ids.size() == 0)
        LOGGER.e("No individuals remaining after applying keep/remove individual files");

    // nothing changed
    if (fam_indi_num == all_ids.size()) {
        valid_indi_num = fam_indi_num;
        genotype.initialize_mask(fam_indi_num);
        LOGGER.i("individuals will be used for analysis", to_string(fam_indi_num));
        return;
    }

    // actually all_ids is already sorted, but just to be safe
    sort(all_ids.begin(), all_ids.end());
    valid_indi_num = all_ids.size();
    
    // get binary mask for genotype reading   
    string FID, IID, str_buf;
    int bit_index = 0, word_index = 0;
    uint64_t word_mask = 0ULL;

    // a coarse way to balance time and memory
    ifstream Fam(famFile.c_str());

    if (valid_indi_num > fam_indi_num / 2) {
        LOGGER.i("Read all individuals into memory and mask unwanted individuals");
        genotype.initialize_mask(fam_indi_num);
        
        while (getline(Fam, str_buf)) {
            if (str_buf.empty()) continue;

            const char* pt = str_buf.c_str();
            parse_string(pt, FID);
            parse_string(pt, IID);

            if (binary_search(all_ids.begin(), all_ids.end(), FID+':'+IID))
                word_mask |= (1ULL << bit_index);

            bit_index++;
            if (bit_index == 64) {
                genotype.X_mask[word_index] = word_mask;
                word_mask = 0ULL;
                bit_index = 0;
                word_index++;
            }
        }

        if (bit_index > 0)
            genotype.X_mask[word_index] = word_mask;
    } else {
        LOGGER.i("Only read individuals to be used into memory");
        genotype.initialize_mask(valid_indi_num);
        
        while (getline(Fam, str_buf)) {
            if (str_buf.empty()) continue;

            const char* pt = str_buf.c_str();
            parse_string(pt, FID);
            parse_string(pt, IID);

            if (binary_search(all_ids.begin(), all_ids.end(), FID+':'+IID)) 
                fam_keep_list.push_back(bit_index);

            bit_index++;
        }
    }
    
    Fam.close();
    LOGGER.i("individuals will be used for analysis", to_string(valid_indi_num));
}


void Cohort::read_frq() 
{   
    string frqFile = params.bfile_list[cohort_index]+".frq";
    ifstream Frq(frqFile.c_str());
    if (!Frq) {
        LOGGER.i("Cannot open PLINK FRQ file [" + frqFile + "] to read, skip allele frequency check");
        return;
    }

    LOGGER.i("Reading PLINK FRQ file from [" + frqFile + "] ...");
    string SNP_buf, A1_buf, A2_buf, str_buf;
    double freq;

    vector<string> vs_buf;
    getline(Frq, str_buf); // the header line

    istringstream iss(str_buf);
    string token;
    while (iss >> token) vs_buf.push_back(token);
    if (vs_buf.size() != 6 || vs_buf[0] != "CHR" || vs_buf[1] != "SNP" || 
        vs_buf[2] != "A1" || vs_buf[3] != "A2" || vs_buf[4] != "MAF" || vs_buf[5] != "NCHROBS")
        LOGGER.e("Format error in frq file, please check");

    int ref_index;

    while (getline(Frq, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        parse_string(pt, str_buf);
        parse_string(pt, SNP_buf);

        auto iter = shared.goodSNP_index_map.find(SNP_buf);
        if (iter == shared.goodSNP_index_map.end()) continue;

        ref_index = iter->second;
    
        parse_string(pt, A1_buf, true);
        parse_string(pt, A2_buf, true);
        parse_double(pt, freq);

        if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
        else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {freq = 1 - freq;}
        else {
            LOGGER.w("removed, different A1 and A2 from ref BIM file, please check", SNP_buf);
            shared.goodSNP_index_map.erase(iter);
            continue;
        } 

        if (abs(sumstat(ref_index, 3) - freq) > params.diff_freq) {
            // LOGGER.w("removed, allele frequency too different between sumstat and bedfile", SNP_buf);
            shared.goodSNP_index_map.erase(SNP_buf);
            continue;
        }
    }

    Frq.close();
}


void Cohort::read_bim()
{
    string bimFile = params.bfile_list[cohort_index]+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e("Cannot open PLINK BIM file [" + bimFile + "] to read");
    LOGGER.i("Reading PLINK BIM file from [" + bimFile + "] ...");

    int chr_buf, ibuf, ref_index;
    string SNP_buf, A1_buf, A2_buf, str_buf;

    swap_array.assign(shared.total_SNP_num, false);
    
    while (getline(Bim, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        parse_int(pt, chr_buf);
        parse_string(pt, SNP_buf);

        bim_SNP_list.push_back(SNP_buf);

        auto iter = shared.goodSNP_index_map.find(SNP_buf);
        if (iter == shared.goodSNP_index_map.end()) continue;
        
        ref_index = iter->second;

        // skip cm position
        skip_delim(pt);
        skip_token(pt);

        if (cohort_index == 0) {
            // first cohort, set as reference
            shared.chr_ref[ref_index] = chr_buf;
            parse_int(pt, shared.SNP_pos_ref[ref_index]);
            parse_string(pt, shared.A1_ref[ref_index], true);
            parse_string(pt, shared.A2_ref[ref_index], true);
        } else {
            parse_int(pt, ibuf);
            parse_string(pt, A1_buf, true);
            parse_string(pt, A2_buf, true);

            // other cohorts, check consistency
            if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
            else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {swap_array[ref_index] = true;}
            else {
                LOGGER.w("removed, different A1 and A2 between two BIM files, please check", SNP_buf);
                shared.goodSNP_index_map.erase(iter);
                continue;
            }

            if (chr_buf != shared.chr_ref[ref_index]) {
                LOGGER.w("removed, different chromosome between two BIM files, please check", SNP_buf);
                shared.goodSNP_index_map.erase(iter);
                continue;
            }

            if (ibuf != shared.SNP_pos_ref[ref_index]) {
                LOGGER.w("removed, different SNP position between two BIM files, please check", SNP_buf);
                shared.goodSNP_index_map.erase(iter);
                continue;
            }      
        }
    }

    Bim.close();
}


void Cohort::read_bed() 
{   
    string bedFile = params.bfile_list[cohort_index]+".bed";
    ifstream Bed(bedFile.c_str(), ios::binary | ios::ate);
    if (!Bed) LOGGER.e("Cannot open PLINK BED file [" + bedFile + "] to read");
    LOGGER.i("Reading PLINK BED file from [" + bedFile + "] in SNP-major format ...");

    const int bytes_per_snp = (fam_indi_num + 3) / 4;
    if (Bed.tellg() != streamoff(uint64_t(bytes_per_snp) * bim_SNP_list.size() + 3)) 
        LOGGER.e("PLINK BED file [" + bedFile + "] size does not match .bim and .fam files, please check");

    Bed.seekg(0, ios::beg);
    char ch[3];
    for (int i=0; i<3; i++) {Bed.read(&ch[i], 1);} 
    if (!Bed || ch[0] != 0x6C || ch[1] != 0x1B || ch[2] != 0x01)
        LOGGER.e("PLINK BED file [" + bedFile + "] not in SNP-major mode, please check");

    genotype.resize(shared.total_SNP_num);

    vector<char> buffer(bytes_per_snp);
    
    if (params.thread_num == 1) {
        for (string SNP_buf : bim_SNP_list) {
            Bed.read(buffer.data(), bytes_per_snp);

            auto iter = shared.goodSNP_index_map.find(SNP_buf);
            if (iter == shared.goodSNP_index_map.end()) continue;

            int ref_index = iter->second;
            
            if (fam_keep_list.empty())
                genotype.decode_single_genotype(buffer, ref_index, swap_array[ref_index]);
            else
                genotype.decode_single_genotype(buffer, ref_index, swap_array[ref_index], fam_keep_list);
        }
    } else {
        #pragma omp parallel
        {   
            #pragma omp single nowait
            {   
                for (string SNP_buf : bim_SNP_list) {
                    Bed.read(buffer.data(), bytes_per_snp);

                    auto iter = shared.goodSNP_index_map.find(SNP_buf);
                    if (iter == shared.goodSNP_index_map.end()) continue;

                    int ref_index = iter->second;

                    auto local_buffer = make_shared<vector<char>>(buffer);

                    // spawn decode task
                    #pragma omp task firstprivate(ref_index, local_buffer) 
                    {
                        if (fam_keep_list.empty())
                            genotype.decode_single_genotype(*local_buffer, ref_index, swap_array[ref_index]);
                        else
                            genotype.decode_single_genotype(*local_buffer, ref_index, swap_array[ref_index], fam_keep_list);
                    }
                }

                #pragma omp taskwait
            } 
        }
    }

    double missing_threshold = valid_indi_num * (1 - params.missingness);

    for (auto iter = shared.goodSNP_index_map.begin(); iter != shared.goodSNP_index_map.end(); ) {
        int ref_index = iter->second;

        if (genotype.X_non_NA_indi_num[ref_index] < missing_threshold) {
            // LOGGER.w("removed, missingness too high or all values are identical in bedfile", iter->first);
            iter = shared.goodSNP_index_map.erase(iter);
            continue;
        }
        
        if (abs(sumstat(ref_index, 3) - genotype.X_avg[ref_index]/2) > params.diff_freq) {
            // LOGGER.w("removed, allele frequency too different between sumstat and bedfile", iter->first);
            iter = shared.goodSNP_index_map.erase(iter);
            continue;
        }

        iter++;
    }

    Bed.close();
    LOGGER.i("Finished reading PLINK BED file");
}


void Cohort::read_PLINK_LD() 
{   
    string ldFile = params.bfile_list[cohort_index]+".ld";
    ifstream Ld(ldFile.c_str());
    if (!Ld) LOGGER.e("Cannot open PLINK LD file [" + ldFile + "] to read");
    LOGGER.i("Reading PLINK LD file from [" + ldFile + "] ...");

    LD_matrix.resize(shared.total_SNP_num);

    // Read ld file
    string SNP1_buf, SNP2_buf, r_buf;
    int SNP1_index, SNP2_index;
    
    while (Ld >> SNP1_buf >> SNP2_buf >> r_buf) {
        auto iter1 = shared.goodSNP_index_map.find(SNP1_buf);
        if (iter1 == shared.goodSNP_index_map.end()) continue;
        SNP1_index = iter1->second;

        auto iter2 = shared.goodSNP_index_map.find(SNP2_buf);
        if (iter2 == shared.goodSNP_index_map.end()) continue;
        SNP2_index = iter2->second;

        LD_matrix(SNP1_index, SNP2_index) = atof(r_buf.c_str()); 
    }

    Ld.close();

    if (cohort_index != 0) {
        // adjust LD_packed according to allele coding
        for (int i = 0; i < shared.total_SNP_num; i++) {
            for (int j = i; j < shared.total_SNP_num; j++) {
                if (swap_array[i] ^ swap_array[j]) {
                    LD_matrix(i,j) = -LD_matrix(i,j);
                    LOGGER << "swap " << i << " " << j << endl;
                }
            }
        }
    }

    LOGGER.i("Finished reading PLINK LD file");
}


void MACOJO::read_input_files() 
{   
    vector<string> commonSNP_all_cohorts;

    // if user provides extract SNP file, read file and initialize common SNPs
    if (!params.extract_file.empty()) {
        skim_file(params.extract_file, commonSNP_all_cohorts, 1, false);
        if (commonSNP_all_cohorts.size() == 0)
            LOGGER.e("Input extract SNP file has no SNPs", params.extract_file);

        LOGGER.i("SNPs initialized by the user\n", to_string(commonSNP_all_cohorts.size()));
    }

    // Step 1: get common SNPs across all cohorts
    for (int n = 0; n < cohorts.size(); n++) {
        auto start = chrono::steady_clock::now();

        vector<string> SNP_sumstat, SNP_PLINK, temp;

        string sumstat_file = params.cojo_file_list[n], PLINK_file = params.bfile_list[n];
        skim_file(PLINK_file+".bim", SNP_PLINK, 2, false);
        LOGGER.i("SNPs in [" + PLINK_file + ".bim]", to_string(SNP_PLINK.size()));

        if (n == 0 && commonSNP_all_cohorts.size() == 0)
            commonSNP_all_cohorts = SNP_PLINK;
        else {
            set_intersection(commonSNP_all_cohorts.begin(), commonSNP_all_cohorts.end(), 
                SNP_PLINK.begin(), SNP_PLINK.end(), back_inserter(temp));
            commonSNP_all_cohorts.swap(temp);
        }

        vector<string>().swap(temp); // clear temp
        skim_file(sumstat_file, SNP_sumstat, 1, true);
        LOGGER.i("SNPs in [" + sumstat_file + "]", to_string(SNP_sumstat.size()));

        set_intersection(commonSNP_all_cohorts.begin(), commonSNP_all_cohorts.end(), 
            SNP_sumstat.begin(), SNP_sumstat.end(), back_inserter(temp));
        commonSNP_all_cohorts.swap(temp);
        
        auto end = chrono::steady_clock::now();
        LOGGER.i("common SNPs after including Cohort "+to_string(n+1), to_string(commonSNP_all_cohorts.size()));
        LOGGER << "Time taken: " << chrono::duration<double>(end-start).count() << " seconds" << endl << endl;
    }

    // if user provides exclude SNP file, read file and delete from common SNPs
    if (!params.exclude_file.empty()) {
        vector<string> exclude_SNP, temp;
        skim_file(params.exclude_file, exclude_SNP, 1, false);

        set_difference(commonSNP_all_cohorts.begin(), commonSNP_all_cohorts.end(), 
            exclude_SNP.begin(), exclude_SNP.end(), back_inserter(temp));
        commonSNP_all_cohorts.swap(temp);
        LOGGER.i("common SNPs after excluding user-specified SNPs\n", to_string(commonSNP_all_cohorts.size()));
    }

    if (commonSNP_all_cohorts.size() == 0)
        LOGGER.e("Input data has no common SNPs across all cohorts");
    
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

    // Step 3: read bim file and sumstat files for all cohorts
    // set cohort 1 as reference
    for (auto& c : cohorts) {
        c.read_bim();
        c.read_sumstat();
    }

    // exclude rare SNPs based on freq_threshold
    for (auto iter = shared.goodSNP_index_map.begin(); iter != shared.goodSNP_index_map.end(); ) {
        int ref_index = iter->second;
        bool erase_flag = !params.if_freq_mode_and;

        if (params.if_freq_mode_and) {
            for (auto &c : cohorts)
                erase_flag |= (c.sumstat(ref_index, 3) < params.maf || c.sumstat(ref_index, 3) > 1-params.maf);
        } else {
            for (auto &c : cohorts) 
                erase_flag &= (c.sumstat(ref_index, 3) < params.maf || c.sumstat(ref_index, 3) > 1-params.maf);
        } 
        
        for (auto &c : cohorts)
            erase_flag |= (c.sumstat(ref_index, 6) > 1e10);

        if (erase_flag) {
            // LOGGER.w("removed, rare SNP based on frequency threshold", iter->first);
            iter = shared.goodSNP_index_map.erase(iter);
        }
        else 
            iter++;
    }

    LOGGER.i("common SNPs after removing rare SNPs", to_string(shared.goodSNP_index_map.size()));

    auto end = chrono::steady_clock::now();
    LOGGER << "Time taken: " << chrono::duration<double>(end-start).count() << " seconds" << endl << endl;

    if (shared.goodSNP_index_map.size() == 0)
        LOGGER.e("Input data has no common SNPs");
    
    // Step 4: read bim and bed files for all cohorts
    for (auto& c : cohorts) {
        // bim file of cohort 1 is reference, do not need to check allele
        start = chrono::steady_clock::now();
    
        if (params.if_LD_mode) {
            c.read_frq();
            c.read_PLINK_LD();
        } else {
            c.read_fam();
            c.read_bed();
        }

        end = chrono::steady_clock::now();
        LOGGER.i("common SNPs remaining", to_string(shared.goodSNP_index_map.size()));
        LOGGER << "Time taken: " << chrono::duration<double>(end-start).count() << " seconds" << endl << endl;

        if (shared.goodSNP_index_map.size() == 0)
            LOGGER.e("Input data has no common SNPs");
    }

    // Step 5: finalize SNP index lists
    for (auto iter = shared.goodSNP_index_map.begin(); iter != shared.goodSNP_index_map.end(); iter++) 
        screened_SNP.push_back(iter->second);

    vector<int> count_array(shared.total_SNP_num);
    for (int i = 0; i < shared.total_SNP_num; i++)
        count_array[i] = i;
    
    set_difference(count_array.begin(), count_array.end(), 
        screened_SNP.begin(), screened_SNP.end(), back_inserter(bad_SNP));

    LOGGER.i("common SNPs at last for analysis", to_string(screened_SNP.size()));
    LOGGER << "--------------------------------" << endl << endl;
}
