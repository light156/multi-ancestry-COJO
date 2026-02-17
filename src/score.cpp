#include "include/macojo.h"


// Compute polygenic scores in a PLINK-compatible --score mode.
void Cohort::calc_polygenic_score(int argc, char** argv) 
{    
    // Step 0: parse command line options
    CLI::App app;

    app.get_formatter()->column_width(50);
    app.get_formatter()->right_column_width(75);

    vector<string> score_options;
    vector<string> extract_options, exclude_options, extract_SNPs, exclude_SNPs;
    string bfile, keep_file, remove_file, output_name;
    int thread_num, bed_block_mb;

    // Main logic options
    app.add_option("--score", score_options, "Per-SNP score file with at least three columns")->required();
    app.add_option("--bfile", bfile, "PLINK binary file prefix [.bim .bed .fam]")->required();
    app.add_option("--out", output_name, "Output file path")->required();
    app.add_option("--keep", keep_file, "File path of individuals to be included");
    app.add_option("--remove", remove_file, "File path of individuals to be excluded");
    app.add_option("--extract", extract_options, "File path of SNPs to be included");
    app.add_option("--exclude", exclude_options, "File path of SNPs to be excluded");
    app.add_option("--extract-snp", extract_SNPs, "A list of SNPs to be included");
    app.add_option("--exclude-snp", exclude_SNPs, "A list of SNPs to be excluded");
    app.add_option("--thread-num", thread_num, "Number of threads to use")->default_val(1)->check(CLI::PositiveNumber);
    app.add_option("--bed-block-mb", bed_block_mb, "BED read block size in MB")->default_val(64)->check(CLI::PositiveNumber);
    
    for (auto *opt : app.get_options())
        opt->multi_option_policy(CLI::MultiOptionPolicy::Throw);

    try { 
        app.parse(argc, argv); 
    } catch (const CLI::CallForHelp &e) {
        LOGGER.i(app.help());
        return;
    } catch (const CLI::ParseError &e) {
        LOGGER.i(app.help());
        LOGGER.e(e.what());
    }
    
    LOGGER.open(output_name + ".prs.log");

    params.bfile_list.push_back(bfile);
    params.keep_file_list.push_back(keep_file);
    params.remove_file_list.push_back(remove_file);
    
    // set number of threads
    #if HAS_OPENMP 
        omp_set_num_threads(thread_num);
    #else
        if (thread_num > 1) {
            LOGGER.w("OpenMP is not available, running in single-thread mode");
            thread_num = 1;
        }
    #endif

    if (!extract_options.empty() || !extract_SNPs.empty()) {
        skim_SNP(extract_options, extract_SNPs);
        if (extract_SNPs.size() == 0)
            LOGGER.e("Please provide valid SNPs for --extract and --extract-snp");

        LOGGER.i("user-specified SNPs to include", extract_SNPs.size());
    }

    if (!exclude_options.empty() || !exclude_SNPs.empty()) {
        skim_SNP(exclude_options, exclude_SNPs);
        if (exclude_SNPs.size() == 0)
            LOGGER.e("Please provide valid SNPs for --exclude and --exclude-snp");
            
        LOGGER.i("user-specified SNPs to exclude", exclude_SNPs.size());
    }

    // Step 1: read score file
    auto start = steady_clock::now();

    string filename = score_options[0];
    vector<string> score_SNP, score_allele, bim_SNP;
    vector<double> score_list;

    int SNP_col = 1;
    int allele_col = 2;
    int score_col = 3;
    bool has_header = false;
    bool output_sum = false;

    vector<string> options_copy = score_options;
    options_copy.erase(options_copy.begin()); // remove filename

    auto before_size = options_copy.size();
    options_copy.erase(remove(options_copy.begin(), options_copy.end(), "header"), options_copy.end());
    if (options_copy.size() < before_size) has_header = true;

    before_size = options_copy.size();
    options_copy.erase(remove(options_copy.begin(), options_copy.end(), "sum"), options_copy.end());
    if (options_copy.size() < before_size) output_sum = true;

    if (!options_copy.empty() && options_copy.size() != 3)
        LOGGER.e("Invalid format for providing score file");

    if (options_copy.size() == 3) {
        const char* pt_SNP = options_copy[0].c_str();
        const char* pt_allele = options_copy[1].c_str();
        const char* pt_score = options_copy[2].c_str();
        if (!parse_num(pt_SNP, SNP_col) || SNP_col < 1)
            LOGGER.e("Invalid SNP column number, should be an integer starting from 1", filename);
        if (!parse_num(pt_allele, allele_col) || allele_col < 1)
            LOGGER.e("Invalid A1 column number, should be an integer starting from 1", filename);
        if (!parse_num(pt_score, score_col) || score_col < 1)
            LOGGER.e("Invalid bJ column number, should be an integer starting from 1", filename);
    }

    if (SNP_col == allele_col || SNP_col == score_col || allele_col == score_col)
        LOGGER.e("SNP, A1, and bJ column numbers must be distinct", filename);

    SNP_col -= 1; // convert to 0-based index
    allele_col -= 1;
    score_col -= 1;
    const int max_col = max(SNP_col, max(allele_col, score_col));

    ifstream sFile(filename.c_str());
    if (!sFile) LOGGER.e("Cannot open [" + filename + "] to read");

    string str_buf;
    if (has_header) getline(sFile, str_buf); // skip header line

    while (getline(sFile, str_buf)) {
        if (str_buf.empty()) continue;

        const char* pt = str_buf.c_str();
        string SNP_buf, allele_buf;
        double score_val = 0.0;
        int ref_index;

        for (int col = 0; col <= max_col; col++) {
            skip_delim(pt);
            if (!*pt) LOGGER.e("Requested column exceeds the total number of columns", filename);

            const char* start = pt;
            skip_token(pt);

            if (col == SNP_col) {
                SNP_buf.assign(start, pt - start);
            } else if (col == allele_col) {
                allele_buf.assign(start, pt - start);
                for (char& c : allele_buf) {
                    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 32);
                }
            } else if (col == score_col) {
                const char* num_pt = start;
                if (!parse_num(num_pt, score_val))
                    LOGGER.e("Invalid numeric value in score file", filename);
            }
        }

        score_SNP.push_back(SNP_buf);
        score_allele.push_back(allele_buf);
        score_list.push_back(score_val);
    }

    sFile.close();

    vector<string> common_SNP(score_SNP);
    sort(common_SNP.begin(), common_SNP.end());
    common_SNP.erase(unique(common_SNP.begin(), common_SNP.end()), common_SNP.end());

    if (common_SNP.size() != score_SNP.size())
        LOGGER.w("Duplicate SNPs found in score file, please make sure this is intended");

    LOGGER.i("Total SNPs in score file", common_SNP.size());

    // user provides extract SNPs
    if (!extract_SNPs.empty()) {
        vector<string> temp;
        set_intersection(common_SNP.begin(), common_SNP.end(), extract_SNPs.begin(), extract_SNPs.end(), back_inserter(temp));
        common_SNP.swap(temp);
        LOGGER.i("common SNPs after including user-specified SNPs\n", common_SNP.size());
    }

    // user provides exclude SNPs
    if (!exclude_SNPs.empty()) {
        vector<string> temp;
        set_difference(common_SNP.begin(), common_SNP.end(), exclude_SNPs.begin(), exclude_SNPs.end(), back_inserter(temp));
        common_SNP.swap(temp);
        LOGGER.i("common SNPs after excluding user-specified SNPs\n", common_SNP.size());
    }
    
    if (common_SNP.size() == 0) LOGGER.e("No SNPs remaining in score file after applying user-specified filters");

    vector<pair<string, int>> score_SNP_table;
    score_SNP_table.reserve(common_SNP.size());

    for (int i = 0; i < score_SNP.size(); i++) {
        if (binary_search(common_SNP.begin(), common_SNP.end(), score_SNP[i]))
            score_SNP_table.emplace_back(score_SNP[i], i);
    }

    sort(score_SNP_table.begin(), score_SNP_table.end(),
        [](const pair<string, int>& a, const pair<string, int>& b) { return a.first < b.first; });

    // Step 2: read bim file and match allele
    vector<int> good_row;
    vector<bool> bed_swap;
    vector<double> good_row_score;

    string bimFile = bfile + ".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e("Cannot open PLINK BIM file [" + bimFile + "] to read");
    LOGGER.i("Reading PLINK BIM file from [" + bimFile + "] ...");

    string SNP_buf, A1_buf, A2_buf;
    int line_num = 0;

    while (getline(Bim, str_buf)) {
        if (str_buf.empty()) continue;
        line_num++;

        const char* pt = str_buf.c_str();
        // skip chr
        skip_delim(pt);
        skip_token(pt);
        parse_string(pt, SNP_buf);

        auto iter = fast_lookup(score_SNP_table, SNP_buf);
        if (iter == score_SNP_table.end()) continue;
        int ref_index = iter->second;

        // skip cm position
        skip_delim(pt);
        skip_token(pt);
        // skip bp
        skip_delim(pt);
        skip_token(pt);
        parse_string(pt, A1_buf, true);
        parse_string(pt, A2_buf, true);

        if (A1_buf == score_allele[ref_index]) {
            good_row.push_back(line_num - 1);
            bed_swap.push_back(false);
            good_row_score.push_back(score_list[ref_index]);
        } else if (A2_buf == score_allele[ref_index]) {
            good_row.push_back(line_num - 1);
            bed_swap.push_back(true);
            good_row_score.push_back(score_list[ref_index]);
        } else {
            LOGGER.w("Allele mismatch for SNP", SNP_buf);
            continue;
        }
    }

    Bim.close();

    LOGGER.i("common SNPs after reading BIM file", good_row.size());
    if (good_row.size() == 0) LOGGER.e("No common SNPs found after reading BIM file");

    // get participants in fam file
    read_fam();
        
    auto end = steady_clock::now();
    LOGGER << "Time taken: " << duration<double>(end-start).count() << " seconds" << endl << endl;

    // Step 3: read bed file and calculate PRS
    start = steady_clock::now();

    string bedFile = bfile + ".bed";
    ifstream Bed(bedFile.c_str(), ios::binary | ios::ate);
    if (!Bed) LOGGER.e("Cannot open PLINK BED file [" + bedFile + "] to read");
    LOGGER.i("Reading PLINK BED file from [" + bedFile + "] in SNP-major format ...");

    const int bytes_per_snp = (fam_indi_num + 3) / 4;
    if (Bed.tellg() != streamoff(uint64_t(bytes_per_snp) * line_num + 3)) 
        LOGGER.e("PLINK BED file [" + bedFile + "] size does not match .bim and .fam files, please check");

    Bed.seekg(0, ios::beg);
    char ch[3];
    for (int i=0; i<3; i++) {Bed.read(&ch[i], 1);} 
    if (!Bed || ch[0] != 0x6C || ch[1] != 0x1B || ch[2] != 0x01)
        LOGGER.e("PLINK BED file [" + bedFile + "] not in SNP-major mode, please check");

    vector<vector<double>> geno_tls(thread_num, vector<double>(fam_indi_num, 0.0));
    vector<vector<int>> cnt_tls(thread_num, vector<int>(fam_indi_num, 0));
    vector<vector<int>> cnt2_tls(thread_num, vector<int>(fam_indi_num, 0));

    const size_t block_bytes = size_t(bed_block_mb) * 1024ULL * 1024ULL;
    const size_t block_snps = std::max<size_t>(1, block_bytes / bytes_per_snp);
    vector<char> block_buffer(block_snps * bytes_per_snp);

    size_t j = 0;
    while (j < good_row.size()) {
        int block_start = good_row[j];
        size_t max_snps = std::min(block_snps, size_t(line_num - block_start));
        size_t bytes_to_read = max_snps * bytes_per_snp;

        Bed.seekg(3 + uint64_t(block_start) * bytes_per_snp, ios::beg);
        Bed.read(block_buffer.data(), bytes_to_read);
        if (!Bed || Bed.gcount() != static_cast<std::streamsize>(bytes_to_read))
            LOGGER.e("Failed to read expected BED block data");

        int block_end = block_start + static_cast<int>(max_snps);
        size_t start_idx = j;
        while (j < good_row.size() && good_row[j] < block_end)
            j++;
        size_t end_idx = j;

        #pragma omp parallel for schedule(static)
        for (size_t t = start_idx; t < end_idx; t++) {
            size_t offset = size_t(good_row[t] - block_start) * bytes_per_snp;
            size_t thread_id = omp_get_thread_num();
            genotype.calc_single_genotype_prs(block_buffer.data() + offset, bytes_per_snp,
                bed_swap[t], good_row_score[t], geno_tls[thread_id], cnt_tls[thread_id], cnt2_tls[thread_id]);
        }
    }

    Bed.close();

    VectorXd geno_vec = VectorXd::Zero(fam_indi_num);
    VectorXi geno_count_vec = VectorXi::Zero(fam_indi_num);
    VectorXi geno_count2_vec = VectorXi::Zero(fam_indi_num);

    for (int t = 0; t < thread_num; t++) {
        geno_vec += VectorXd::Map(geno_tls[t].data(), fam_indi_num);
        geno_count_vec += VectorXi::Map(cnt_tls[t].data(), fam_indi_num);
        geno_count2_vec += VectorXi::Map(cnt2_tls[t].data(), fam_indi_num);
    }

    ofstream PRS_file((output_name + ".prs.profile").c_str());
    PRS_file << "FID\tIID\tCNT\tCNT2\tSCORE" << endl;

    int total_allele_num = good_row.size() * 2;

    for (size_t w = 0; w < genotype.words_per_snp; w++) {
        uint64_t mask = genotype.X_mask[w];
        while (mask != 0ULL) {
            int bit = __builtin_ctzll(mask);
            int n = static_cast<int>(w * 64 + bit);
            string s = fam_ID_array[n];
            auto pos = s.find(':');
            double score_out = output_sum ? geno_vec(n) : (geno_vec(n) / total_allele_num);
            PRS_file << s.substr(0, pos) << "\t" << s.substr(pos + 1) << "\t" 
                << total_allele_num + geno_count_vec(n) << "\t" 
                << geno_count2_vec(n) << "\t" 
                << score_out << endl;
            mask &= (mask - 1ULL);
        }
    }

    PRS_file.close();
    LOGGER.i("Polygenic scores written to file [" + output_name + ".prs.profile]");

    end = steady_clock::now();
    LOGGER << "Time taken: " << duration<double>(end-start).count() << " seconds" << endl << endl;
}
