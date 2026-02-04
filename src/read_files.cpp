#include "include/macojo.h"


// Parse SNP list files or inline SNP names with optional header/column index.
void skim_SNP(const vector<string>& options, vector<string>& SNP_list)
{   
    if (!options.empty()) 
    {
        string filename = options[0];
        int col_idx = 1;
        bool has_header = false;

        vector<string> options_copy = options;
        options_copy.erase(options_copy.begin()); // remove filename
        
        options_copy.erase(remove(options_copy.begin(), options_copy.end(), "header"), options_copy.end());
        if (options_copy.size() < options.size()) has_header = true;

        if (options_copy.size() > 1) LOGGER.e("Invalid format for providing SNP file");
        
        if (!options_copy.empty()) {
            const char* pt = options_copy[0].c_str();
            if (!parse_num(pt, col_idx) || col_idx < 1) 
                LOGGER.e("Invalid column number, should be an integer starting from 1", filename);
        } 

        col_idx -= 1; // convert to 0-based index

        ifstream sFile(filename.c_str());
        if (!sFile) LOGGER.e("Cannot open [" + filename + "] to read");

        string str_buf;
        if (has_header) getline(sFile, str_buf); // skip header line

        while (getline(sFile, str_buf)) {
            if (str_buf.empty()) continue;

            const char* pt = str_buf.c_str();
            
            // iterate to column col_idx
            int col = 0;
            while (col < col_idx) {
                skip_delim(pt);
                skip_token(pt);
                if (!*pt) LOGGER.e("Requested column exceeds the total number of columns", filename);
                col++;
            }

            skip_delim(pt);

            // p now at start of desired column
            const char* start = pt;
            skip_token(pt);
            SNP_list.emplace_back(start, pt - start);
        }

        sFile.close();
    }

    sort(SNP_list.begin(), SNP_list.end());

    for (int i = 1; i < SNP_list.size(); i++) {
        if (SNP_list[i] == SNP_list[i-1]) {
            LOGGER.w("Duplicate SNP found", SNP_list[i]);
            while (i+1 < SNP_list.size() && SNP_list[i+1] == SNP_list[i]) i++; // skip long runs of the same duplicate
        }
    }

    SNP_list.erase(unique(SNP_list.begin(), SNP_list.end()), SNP_list.end());
    SNP_list.erase(remove(SNP_list.begin(), SNP_list.end(), ""), SNP_list.end());
    SNP_list.erase(remove(SNP_list.begin(), SNP_list.end(), "."), SNP_list.end());
}


// Read FAM file and emit "FID:IID" identifiers.
void skim_fam(string filename, vector<string>& str_list) 
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
}


// Read BIM file and extract SNP names for a given chromosome.
void skim_bim(string filename, int chr, vector<string>& SNP_list)
{
    ifstream sFile(filename.c_str());
    if (!sFile) LOGGER.e("Cannot open file [" + filename + "] to read");

    string SNP_buf, str_buf;
    int chr_buf;

    while (getline(sFile, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        parse_num(pt, chr_buf);
        if (chr_buf != chr) continue;

        parse_string(pt, SNP_buf);
        SNP_list.push_back(SNP_buf);
    }
    
    sFile.close();
    
    sort(SNP_list.begin(), SNP_list.end());

    for (int i = 1; i < SNP_list.size(); i++) {
        if (SNP_list[i] == SNP_list[i-1]) {
            LOGGER.w("Duplicate SNP found", SNP_list[i]);
            while (i+1 < SNP_list.size() && SNP_list[i+1] == SNP_list[i]) i++; // skip long runs of the same duplicate
        }
    }

    SNP_list.erase(unique(SNP_list.begin(), SNP_list.end()), SNP_list.end());
    SNP_list.erase(remove(SNP_list.begin(), SNP_list.end(), ""), SNP_list.end());
    SNP_list.erase(remove(SNP_list.begin(), SNP_list.end(), "."), SNP_list.end());
}


// Read FAM and apply keep/remove filters to build a mask.
void Cohort::read_fam() 
{   
    string famFile = params.bfile_list[cohort_index]+".fam";
    skim_fam(famFile, fam_ID_array);
    
    fam_indi_num = fam_ID_array.size();
    vector<string> sorted_ID(fam_ID_array);
    sort(sorted_ID.begin(), sorted_ID.end());
    if (adjacent_find(sorted_ID.begin(), sorted_ID.end()) != sorted_ID.end())
        LOGGER.e("Duplicate individuals in [" + famFile + "], please check", *adjacent_find(sorted_ID.begin(), sorted_ID.end()));

    LOGGER.i("individuals in FAM file [" + famFile + "]", fam_indi_num);

    if (!params.keep_file_list[cohort_index].empty()) {
        vector<string> keep_ids, temp_ids;
        skim_fam(params.keep_file_list[cohort_index], keep_ids);        
        
        sort(keep_ids.begin(), keep_ids.end());
        set_intersection(sorted_ID.begin(), sorted_ID.end(), keep_ids.begin(), keep_ids.end(), back_inserter(temp_ids));
        sorted_ID.swap(temp_ids);
        LOGGER.i("individuals after keeping individuals", sorted_ID.size());
    }

    if (!params.remove_file_list[cohort_index].empty()) {
        vector<string> remove_ids, temp_ids;
        skim_fam(params.remove_file_list[cohort_index], remove_ids);

        sort(remove_ids.begin(), remove_ids.end());
        set_difference(sorted_ID.begin(), sorted_ID.end(), remove_ids.begin(), remove_ids.end(), back_inserter(temp_ids));
        sorted_ID.swap(temp_ids);
        LOGGER.i("individuals after removing individuals", sorted_ID.size());
    }

    if (sorted_ID.size() == 0)
        LOGGER.e("No individuals remaining after applying keep/remove individual files");
    
    valid_indi_num = sorted_ID.size();

    // nothing changed
    if (fam_indi_num == valid_indi_num) {
        genotype.initialize_mask(fam_indi_num);
        LOGGER.i("individuals will be used for analysis", fam_indi_num);
        return;
    }

    // get binary mask for genotype reading 
    int bit_index = 0, word_index = 0;
    uint64_t word_mask = 0ULL;

    // a coarse way to balance time and memory
    genotype.initialize_mask(fam_indi_num);
    
    for (const auto& full_ID : fam_ID_array) {
        if (binary_search(sorted_ID.begin(), sorted_ID.end(), full_ID))
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
    
    LOGGER.i("individuals will be used for analysis", valid_indi_num);
}


// Read GWAS summary statistics in GCTA-COJO format.
void Cohort::read_sumstat() 
{  
    string cojo_file = params.cojo_file_list[cohort_index];
    ifstream Meta(cojo_file.c_str());
    if (!Meta) LOGGER.e("Cannot open sumstat file [" + cojo_file + "] to read");
    LOGGER.i("Reading GWAS summary-level statistics from [" + cojo_file + "] ...");

    string SNP_buf, A1_buf, A2_buf, str_buf;
    double freq, b, se, p, N, V;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line

    istringstream iss(str_buf);
    string token;
    while (iss >> token) vs_buf.push_back(token);
    if (vs_buf.size() < 8)
        LOGGER.e("Sumstat file should be in GCTA-COJO format, with a header like: SNP A1 A2 freq b se p N,"
                 "and the first 8 columns should be in this order");
    
    int ref_index;
    vector<double> Vp_gcta_list;

    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D
    sumstat.resize(shared.goodSNP_table.size(), 7);
    
    while (getline(Meta, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        parse_string(pt, SNP_buf);
        parse_string(pt, A1_buf, true);
        parse_string(pt, A2_buf, true);

        auto iter = fast_lookup(shared.goodSNP_table, SNP_buf);

        if (!parse_num(pt, freq)) { 
            if (iter != shared.goodSNP_table.end()) {
                shared.bad_SNP_dict.emplace_back(BadSnpReason::InvalidNumericValue, SNP_buf);
                iter->second = -1;
            }
            continue;
        } 

        if (freq < 1e-10 || freq > 1.0 - 1e-10) { 
            if (iter != shared.goodSNP_table.end()) iter->second = -1;
            Vp_gcta_list.push_back(0);
            continue;
        } 

        if (!parse_num(pt, b) || !parse_num(pt, se) || !parse_num(pt, p) || !parse_num(pt, N)) {
            if (iter != shared.goodSNP_table.end()) {
                shared.bad_SNP_dict.emplace_back(BadSnpReason::InvalidNumericValue, SNP_buf);
                iter->second = -1;
            }
            continue;
        } 

        V = 2 * freq * (1 - freq);
        Vp_gcta_list.push_back(V * N * (se * se + b * b / (N - 1.0)));

        if (iter == shared.goodSNP_table.end()) continue;
        ref_index = iter->second;

        if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
        else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {freq = 1 - freq; b = -b;}
        else {
            shared.bad_SNP_dict.emplace_back(BadSnpReason::AlleleMismatch, SNP_buf);
            iter->second = -1;
            continue;
        } 

        // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D
        sumstat.row(ref_index) << b, se * se, p, freq, N, V, 0.0;
    }

    Meta.close();

    Vp = median(Vp_gcta_list);
    sumstat.col(4) = (Vp / sumstat.col(5) - square(sumstat.col(0))) / sumstat.col(1) + 1;
    sumstat.col(6) = sumstat.col(4) * sumstat.col(5);

    if (params.if_infoscore)
        sumstat.col(0) = sumstat.col(0) * sqrt(sumstat.col(4) / sumstat.col(4).maxCoeff());

    LOGGER << "Estimated phenotype variance Vp: " << Vp << endl;
}


// Optional frequency check against PLINK .frq if present.
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
    int chr_buf, ref_index;
    double freq;

    vector<string> vs_buf;
    getline(Frq, str_buf); // the header line

    istringstream iss(str_buf);
    string token;
    while (iss >> token) vs_buf.push_back(token);
    if (vs_buf.size() < 6 || vs_buf[0] != "CHR" || vs_buf[1] != "SNP" || 
        vs_buf[2] != "A1" || vs_buf[3] != "A2" || vs_buf[4] != "MAF" || vs_buf[5] != "NCHROBS")
        LOGGER.e("Format error in frq file, please check");

    while (getline(Frq, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();

        parse_num(pt, chr_buf);
        if (chr_buf != params.curr_chr) continue;

        parse_string(pt, SNP_buf);
        auto iter = fast_lookup(shared.goodSNP_table, SNP_buf);
        if (iter == shared.goodSNP_table.end()) continue;
        ref_index = iter->second;

        parse_string(pt, A1_buf, true);
        parse_string(pt, A2_buf, true);
        parse_num(pt, freq);

        if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
        else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {freq = 1 - freq;}
        else {
            shared.bad_SNP_dict.emplace_back(BadSnpReason::AlleleMismatch, SNP_buf);
            iter->second = -1;
            continue;
        } 

        if (abs(sumstat(ref_index, 3) - freq) > params.diff_freq) {
            shared.bad_SNP_dict.emplace_back(BadSnpReason::FreqDiffTooLarge, SNP_buf);
            iter->second = -1;
            continue;
        }
    }

    Frq.close();
}


// Read BIM for allele checks and reference positions.
void Cohort::read_bim()
{
    string bimFile = params.bfile_list[cohort_index]+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e("Cannot open PLINK BIM file [" + bimFile + "] to read");
    LOGGER.i("Reading PLINK BIM file from [" + bimFile + "] ...");

    int chr_buf, ibuf, ref_index;
    string SNP_buf, A1_buf, A2_buf, str_buf;

    bed_swap_array.assign(shared.goodSNP_table.size(), false);
    
    while (getline(Bim, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        parse_num(pt, chr_buf);
        parse_string(pt, SNP_buf);
        bim_SNP_array.push_back(SNP_buf);

        if (chr_buf != params.curr_chr) continue;

        auto iter = fast_lookup(shared.goodSNP_table, SNP_buf);
        if (iter == shared.goodSNP_table.end()) continue;
        ref_index = iter->second;
    
        // skip cm position
        skip_delim(pt);
        skip_token(pt);
        parse_num(pt, ibuf);
        parse_string(pt, A1_buf, true);
        parse_string(pt, A2_buf, true);

        // first cohort, set as reference
        if (cohort_index == 0) { 
            shared.SNP_pos_ref[ref_index] = ibuf;
            shared.A1_ref[ref_index] = A1_buf;
            shared.A2_ref[ref_index] = A2_buf;
        } else {
            if (shared.A1_ref[ref_index] == A1_buf && shared.A2_ref[ref_index] == A2_buf) {}
            else if (shared.A1_ref[ref_index] == A2_buf && shared.A2_ref[ref_index] == A1_buf) {bed_swap_array[ref_index] = true;}
            else {
                shared.bad_SNP_dict.emplace_back(BadSnpReason::AlleleMismatch, SNP_buf);
                iter->second = -1;
                continue;
            }

            if (ibuf != shared.SNP_pos_ref[ref_index]) {
                shared.bad_SNP_dict.emplace_back(BadSnpReason::BpMismatch, SNP_buf);
                iter->second = -1;
                continue;
            }      
        }
    }

    Bim.close();
}


// Read PLINK .bed in SNP-major format and populate Geno buffers.
void Cohort::read_bed() 
{   
    string bedFile = params.bfile_list[cohort_index]+".bed";
    ifstream Bed(bedFile.c_str(), ios::binary | ios::ate);
    if (!Bed) LOGGER.e("Cannot open PLINK BED file [" + bedFile + "] to read");
    LOGGER.i("Reading PLINK BED file from [" + bedFile + "] in SNP-major format ...");

    const int bytes_per_snp = (fam_indi_num + 3) / 4;
    if (Bed.tellg() != streamoff(uint64_t(bytes_per_snp) * bim_SNP_array.size() + 3)) 
        LOGGER.e("PLINK BED file [" + bedFile + "] size does not match .bim and .fam files, please check");

    Bed.seekg(0, ios::beg);
    char ch[3];
    for (int i=0; i<3; i++) {Bed.read(&ch[i], 1);} 
    if (!Bed || ch[0] != 0x6C || ch[1] != 0x1B || ch[2] != 0x01)
        LOGGER.e("PLINK BED file [" + bedFile + "] not in SNP-major mode, please check");

    genotype.resize(shared.goodSNP_table.size());

    vector<pair<int, int>> selected;
    for (int j = 0; j < bim_SNP_array.size(); j++) {
        auto iter = fast_lookup(shared.goodSNP_table, bim_SNP_array[j]);
        if (iter == shared.goodSNP_table.end()) continue;
        
        selected.emplace_back(j, iter->second);
    }

    const size_t block_bytes = size_t(params.bed_block_mb) * 1024ULL * 1024ULL;
    const size_t block_snps = std::max<size_t>(1, block_bytes / bytes_per_snp);
    vector<char> block_buffer(block_snps * bytes_per_snp);

    size_t idx = 0;
    while (idx < selected.size()) {
        int block_start = selected[idx].first;
        size_t max_snps = std::min(block_snps, size_t(bim_SNP_array.size() - block_start));
        size_t bytes_to_read = max_snps * bytes_per_snp;

        Bed.seekg(3 + uint64_t(block_start) * bytes_per_snp, ios::beg);
        Bed.read(block_buffer.data(), bytes_to_read);
        if (!Bed || Bed.gcount() != static_cast<std::streamsize>(bytes_to_read))
            LOGGER.e("Failed to read expected BED block data");

        int block_end = block_start + static_cast<int>(max_snps);
        size_t start_idx = idx;
        while (idx < selected.size() && selected[idx].first < block_end)
            idx++;
        size_t end_idx = idx;

        #pragma omp parallel for schedule(static)
        for (size_t t = start_idx; t < end_idx; t++) {
            size_t offset = size_t(selected[t].first - block_start) * bytes_per_snp;
            int ref_index = selected[t].second;
            genotype.decode_single_genotype(block_buffer.data() + offset, bytes_per_snp,
                ref_index, bed_swap_array[ref_index]);
        }
    }

    double missing_threshold = valid_indi_num * (1 - params.missingness);

    for (auto iter = shared.goodSNP_table.begin(); iter != shared.goodSNP_table.end(); iter++) {
        if (iter->second == -1) continue;
        int ref_index = iter->second;

        if (genotype.X_non_NA_indi_num[ref_index] < missing_threshold) {
            shared.bad_SNP_dict.emplace_back(BadSnpReason::HighMissingness, iter->first);
            iter->second = -1;
            continue;
        }
        
        if (abs(sumstat(ref_index, 3) - genotype.X_avg[ref_index]/2) > params.diff_freq) {
            shared.bad_SNP_dict.emplace_back(BadSnpReason::FreqDiffTooLarge, iter->first);
            iter->second = -1;
            continue;
        }
    }

    Bed.close();
    vector<string>().swap(bim_SNP_array);
}


// Read PLINK .ld file into compact upper-triangle storage.
void Cohort::read_PLINK_LD() 
{   
    string ldFile = params.bfile_list[cohort_index]+".ld";
    ifstream Ld(ldFile.c_str());
    if (!Ld) LOGGER.e("Cannot open PLINK LD file [" + ldFile + "] to read");
    LOGGER.i("Reading PLINK LD file from [" + ldFile + "] ...");

    int total_SNP_num = shared.goodSNP_table.size();
    LD_matrix.resize(total_SNP_num);

    // Read ld file
    string SNP1_buf, SNP2_buf, str_buf;
    double r;
    
    while (getline(Ld, str_buf)) {
        if (str_buf.empty()) continue;
        const char* pt = str_buf.c_str();

        parse_string(pt, SNP1_buf);
        auto iter1 = fast_lookup(shared.goodSNP_table, SNP1_buf);
        if (iter1 == shared.goodSNP_table.end()) continue;

        parse_string(pt, SNP2_buf);
        auto iter2 = fast_lookup(shared.goodSNP_table, SNP2_buf);
        if (iter2 == shared.goodSNP_table.end()) continue;

        parse_num(pt, r);
        LD_matrix(iter1->second, iter2->second) = r; 
    }

    Ld.close();

    // bim file of cohort 1 is reference, do not need to check allele
    if (cohort_index != 0) {
        // adjust LD_packed according to allele coding
        for (int i = 0; i < total_SNP_num; i++) {
            for (int j = i; j < total_SNP_num; j++) {
                if (bed_swap_array[i] ^ bed_swap_array[j]) {
                    LD_matrix(i,j) = -LD_matrix(i,j);
                    // LOGGER << "swap " << i << " " << j << endl;
                }
            }
        }
    }

    LOGGER.i("Finished reading PLINK LD file");
}


// Build common SNPs across cohorts, read all inputs, and apply filters.
bool MACOJO::read_input_files() 
{   
    vector<string> common_SNP;

    // Step 1: get common SNPs across all cohorts
    for (int n = 0; n < cohorts.size(); n++) {
        auto start = steady_clock::now();

        vector<string> SNP_PLINK, SNP_sumstat, temp;

        skim_bim(params.bfile_list[n]+".bim", params.curr_chr, SNP_PLINK);    
        LOGGER.i("SNPs in [" + params.bfile_list[n]+".bim" + "]", SNP_PLINK.size());

        if (n == 0)
            common_SNP.swap(SNP_PLINK);
        else {
            set_intersection(common_SNP.begin(), common_SNP.end(), SNP_PLINK.begin(), SNP_PLINK.end(), back_inserter(temp));
            common_SNP.swap(temp);
        }

        skim_SNP({params.cojo_file_list[n], "header"}, SNP_sumstat);
        LOGGER.i("SNPs in [" + params.cojo_file_list[n] + "]", SNP_sumstat.size());

        vector<string>().swap(temp);
        set_intersection(common_SNP.begin(), common_SNP.end(), SNP_sumstat.begin(), SNP_sumstat.end(), back_inserter(temp));
        common_SNP.swap(temp);
        
        auto end = steady_clock::now();
        LOGGER.i("common SNPs after including Cohort " + to_string(n+1), common_SNP.size());
        LOGGER << "Time taken: " << duration<double>(end-start).count() << " seconds" << endl << endl;
    }

    // user provides extract SNPs
    if (!params.extract_SNPs.empty()) {
        vector<string> temp;
        set_intersection(common_SNP.begin(), common_SNP.end(), params.extract_SNPs.begin(), params.extract_SNPs.end(), back_inserter(temp));
        common_SNP.swap(temp);
        LOGGER.i("common SNPs after including user-specified SNPs\n", common_SNP.size());
    }

    // user provides exclude SNPs
    if (!params.exclude_SNPs.empty()) {
        vector<string> temp;
        set_difference(common_SNP.begin(), common_SNP.end(), params.exclude_SNPs.begin(), params.exclude_SNPs.end(), back_inserter(temp));
        common_SNP.swap(temp);
        LOGGER.i("common SNPs after excluding user-specified SNPs\n", common_SNP.size());
    }

    if (common_SNP.size() == 0) return false;
    
    auto start = steady_clock::now();

    int total_SNP_num = common_SNP.size();

    shared.A1_ref.resize(total_SNP_num);
    shared.A2_ref.resize(total_SNP_num);
    shared.SNP_pos_ref.resize(total_SNP_num);
    shared.bp_order.resize(total_SNP_num);
    shared.goodSNP_table.reserve(total_SNP_num);

    for (int i = 0; i < total_SNP_num; i++) {
        shared.bp_order[i] = i;
        shared.goodSNP_table.emplace_back(common_SNP[i], i);
    }

    // Step 2: read bim file and sumstat files for all cohorts, set cohort 1 as reference
    for (auto& c : cohorts) {
        c.read_bim();
        c.read_sumstat();
    }

    // exclude rare SNPs based on freq_threshold
    for (auto iter = shared.goodSNP_table.begin(); iter != shared.goodSNP_table.end(); iter++) {
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
            shared.bad_SNP_dict.emplace_back(BadSnpReason::RareVariant, iter->first);
            iter->second = -1;
        }
    }

    int goodSNP_num = 0;
    for (const auto& kv : shared.goodSNP_table) goodSNP_num += (kv.second != -1);
    LOGGER.i("common SNPs after removing rare SNPs", goodSNP_num);
    if (goodSNP_num == 0) return false;

    auto end = steady_clock::now();
    LOGGER << "Time taken: " << duration<double>(end-start).count() << " seconds" << endl << endl;

    // Step 3: read bed files for all cohorts
    for (auto& c : cohorts) {
        start = steady_clock::now();
    
        if (params.if_LD_mode) {
            c.read_frq();
            c.read_PLINK_LD();
        } else {
            c.read_fam();
            c.read_bed();
        }
        
        end = steady_clock::now();
        LOGGER << "Time taken: " << duration<double>(end-start).count() << " seconds" << endl << endl;
    }

    // finalize SNP index lists; active_mask is initialized once and updated during selection
    active_mask.assign(total_SNP_num, 1);
    for (int i = 0; i < total_SNP_num; i++) {
        if (shared.goodSNP_table[i].second == -1) {
            bad_SNP.push_back(i);
            active_mask[i] = 0;
        }
    }

    if (bad_SNP.size() == total_SNP_num) return false;

    if (params.if_output_all)
        output_bad_SNP(params.output_name);

    // reorder bp_order according to bp position
    sort(shared.bp_order.begin(), shared.bp_order.end(),
        [&](int a, int b) { return shared.SNP_pos_ref[a] < shared.SNP_pos_ref[b]; });

    LOGGER.i("common SNPs at last for analysis", total_SNP_num - bad_SNP.size());
    LOGGER << "--------------------------------" << endl << endl;
    return true;
}
