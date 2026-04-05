#include "include/macojo.h"


// Formatter that suppresses "Excludes: ..." lines from --help output.
// CLI11 appends those automatically when .excludes() is used, cluttering the help.
struct SilentExcludesFormatter : public CLI::Formatter {
    string make_option_opts(const CLI::Option* opt) const override {
        string res = CLI::Formatter::make_option_opts(opt);
        // Strip ":{value1,value2,...}" appended by CLI::IsMember to the type name
        auto brace = res.find(":{");
        if (brace != string::npos) {
            auto close = res.find('}', brace);
            if (close != string::npos) res.erase(brace, close - brace + 1);
        }
        // Strip " Excludes: ..." appended when .excludes() is used
        auto excl = res.find(" Excludes:");
        if (excl != string::npos) res.erase(excl);
        return res;
    }
};


// Parse CLI options, validate combinations, and populate HyperParams.
void set_read_process_output_options(int argc, char** argv)
{   
    HyperParams& params = get_params();
    CLI::App app;

    app.description(
        "\nMulti-ancestry conditional and joint analysis of GWAS summary statistics\n"
        "Code repo: https://github.com/light156/multi-ancestry-COJO\n"
        "Documentation: https://light156.github.io/multi-ancestry-COJO-docs/\n"
        "\nTIPS:\n"
        "  1. The usage is largely consistent with the original GCTA COJO, while extended to multiple cohorts and LD inputs. \n"
        "  2. Please make sure the paths are paired in --bfile/--ld and --cojo-file.\n"
        "     (e.g., --bfile xx1.sumstat xx2.sumstat xx3.sumstat --cojo-file yy1 yy2 yy3)\n"
        "  3. For using genotype data (--bfile), please provide (.bed .bim .fam) files.\n"
        "     For using LD matrix (--ld), please provide at least (.bim .ld) files.\n"
        "  4. The folders in the output file path (which is after --out) must exist.\n"
        "     (You can easily create them with 'mkdir -p <folder_path>')"
    );

    auto fmt = std::make_shared<SilentExcludesFormatter>();
    fmt->column_width(65);
    fmt->right_column_width(75);
    app.formatter(fmt);

    size_t cohort_num;

    // Main logic options
    string io_group = "Input/Output File Options";
    auto *bfile = app.add_option("--bfile", params.bfile_list, "PLINK binary file prefix for each cohort [.bim .bed .fam]")->group(io_group);
    auto *ld =  app.add_option("--ld", params.bfile_list, "PLINK LD file prefix for each cohort [.bim .ld (.frq)]")->group(io_group);
    app.add_option("--cojo-file", params.cojo_file_list, "GWAS summary statistics file for each cohort (GCTA-COJO format)")->required()->group(io_group);
    app.add_option("--out", params.output_name, "Output file path")->required()->group(io_group);

    // options same as GCTA-COJO
    string cojo_group = "Original GCTA Options/Flags";
    auto *slct_flag = app.add_flag("--cojo-slct", "Stepwise iterative selection of independently associated SNPs")->group(cojo_group);
    auto *joint_flag = app.add_flag("--cojo-joint", "Only calculate joint effect for provided SNPs and exit")->group(cojo_group);
    auto *cond_option = app.add_option("--cojo-cond", params.cojo_cond_options, "Only calculate conditional effect for provided SNPs and exit")->group(cojo_group);
    app.add_option("--cojo-wind", params.window_kb, "SNP position window in Kb (-1 for no windows)")->default_val(10000)->group(cojo_group);
    app.add_option("--cojo-p", params.p_value, "Significance threshold for SNP selection")->default_val(5e-8)->check(CLI::Range(0.0, 1.0))->group(cojo_group)->excludes(joint_flag, cond_option);
    app.add_option("--cojo-collinear", params.collinear, "Colinearity threshold for SNP selection")->default_val(0.9)->check(CLI::Range(0.0, 0.9999))->group(cojo_group);
    app.add_option("--cojo-top-SNPs", params.max_iter_num, "Only select a fixed number of independently associated SNPs ")->default_val(10000)->check(CLI::PositiveNumber)->group(cojo_group)->excludes(joint_flag, cond_option);
    app.add_option("--extract", params.extract_options, "File path of SNPs to be included")->group(cojo_group);
    app.add_option("--exclude", params.exclude_options, "File path of SNPs to be excluded")->group(cojo_group);
    app.add_option("--extract-snp", params.extract_SNPs, "A list of SNPs to be included")->group(cojo_group);
    app.add_option("--exclude-snp", params.exclude_SNPs, "A list of SNPs to be excluded")->group(cojo_group);
    app.add_option("--chr", params.curr_chr, "Only include SNPs on a specific chromosome")->group(cojo_group)->check(CLI::Range(1, 22));
    app.add_option("--thread-num", params.thread_num, "Number of threads to use")->default_val(1)->check(CLI::PositiveNumber)->group(cojo_group);
    
    string extend_group = "Multi-ancestry-extended Options";
    app.add_option("--keep", params.keep_file_list, "File paths of individuals to be included")->group(extend_group)->excludes(ld);
    app.add_option("--remove", params.remove_file_list, "File paths of individuals to be excluded")->group(extend_group)->excludes(ld);
    app.add_option("--diff-freq", params.diff_freq_list, "Frequency difference threshold between sumstat and genotype")->check(CLI::Range(0.0, 1.0))->group(extend_group);
    app.add_option("--maf", params.maf_list, "Minor Allele Frequency threshold for genotype")->check(CLI::Range(0.0, 0.5))->group(extend_group);
    app.add_option("--maf-sumstat", params.maf_sumstat_list, "Minor Allele Frequency threshold for sumstat files")->check(CLI::Range(0.0, 0.5))->group(extend_group);
    app.add_option("--geno", params.missingness_list, "Missingness threshold in .bed files (1 for no threshold)")->check(CLI::Range(0.0, 1.0))->group(extend_group)->excludes(ld);
    
    // multi-cohort specific options
    string manc_group = "Manc-COJO Specific Options/Flags";
    app.add_option("--fix", params.fix_options, "File path of fixed SNPs for iterative selection")->group(manc_group)->excludes(joint_flag, cond_option);
    app.add_option("--fix-snp", params.fix_SNPs, "A list of fixed SNPs for iterative selection")->group(manc_group)->excludes(joint_flag, cond_option);
    app.add_option("--fix-p", params.fix_p_value, "Significance threshold for fixed SNPs")->default_val(5e-8)->check(CLI::Range(0.0, 1.0))->group(manc_group)->excludes(joint_flag, cond_option);
    app.add_option("--fix-col", params.fix_collinear, "Collinearity threshold for fixed SNPs")->default_val(0.9)->check(CLI::Range(0.0, 0.9999))->group(manc_group)->excludes(joint_flag, cond_option);
    app.add_option("--fix-cojo-col", params.fix_cojo_collinear, "Collinearity threshold between fixed SNPs and selected SNPs")->default_val(0.9)->check(CLI::Range(0.0, 0.9999))->group(manc_group)->excludes(joint_flag, cond_option);
    app.add_flag("--fix-drop", params.if_fix_drop, "Drop fixed SNPs above --fix-p threshold before SNP selection")->group(manc_group)->excludes(joint_flag, cond_option);
    auto *R2_option     = app.add_option("--R2", params.R2_threshold, "R2 incremental threshold (-1 for no threshold)")->default_val(-1)->group(manc_group)->excludes(joint_flag, cond_option);
    auto *R2back_option = app.add_option("--R2back", params.R2back_threshold, "R2 threshold for backward selection (-1 for no threshold)")->default_val(-1)->group(manc_group)->excludes(joint_flag, cond_option);
    app.add_option("--slct-mode", params.slct_mode, "Iterative SNP selection method (GCTA, imputeNA, or removeNA)")
        ->check(CLI::IsMember({"GCTA", "imputeNA", "removeNA"}))->type_name("MODE")->default_val("GCTA")->group(manc_group)->excludes(joint_flag, cond_option);
    app.add_option("--effect-mode", params.effect_size_mode, "Effect size estimation method (GCTA, imputeNA, or removeNA)")
        ->check(CLI::IsMember({"GCTA", "imputeNA", "removeNA"}))->type_name("MODE")->default_val("GCTA")->group(manc_group);
    app.add_flag("--var-from-ld", params.if_var_from_ld, "Estimate genotypic variance from LD reference instead of 2pq")->group(manc_group)->excludes(ld);
    app.add_flag("--hetero-report", params.if_hetero_report, "Output Cochran's Q heterogeneity statistics (multi-cohort only)")->group(manc_group);
    app.add_flag("--output-all", params.if_output_all, "Save all .cma.cojo .jma.cojo .ldr.cojo and .badsnps files")->group(manc_group);
    app.add_option("--bed-block-mb", params.bed_block_mb, "BED read block size in MB")->default_val(64)->check(CLI::PositiveNumber)->group("");
    
    for (auto *opt : app.get_options())
        opt->multi_option_policy(CLI::MultiOptionPolicy::Throw);

    // Post-parse validation and derived parameter setup.
    app.callback([&]() {
        cohort_num = params.cojo_file_list.size();

        if ((slct_flag->count() > 0) + (joint_flag->count() > 0) + (cond_option->count() > 0) != 1)
            throw CLI::ValidationError("Exactly one of --cojo-slct, --cojo-joint, or --cojo-cond must be provided to specify the analysis mode");
        
        if ((bfile->count() > 0) == (ld->count() > 0))
            throw CLI::ValidationError("Either --bfile or --ld must be provided for all cohorts");

        if (params.bfile_list.size() != cohort_num)
            throw CLI::ValidationError("Filepath numbers after --bfile/--ld and --cojo-file must be the same");

        if (params.keep_file_list.empty())
            params.keep_file_list.assign(cohort_num, "");
        else if (params.keep_file_list.size() != cohort_num)
            throw CLI::ValidationError("Filepath numbers after --cojo-file and --keep must be the same");

        if (params.remove_file_list.empty())
            params.remove_file_list.assign(cohort_num, "");
        else if (params.remove_file_list.size() != cohort_num)
            throw CLI::ValidationError("Filepath numbers after --cojo-file and --remove must be the same");

        auto expand_or_validate = [&](vector<double>& v, double default_val, const string& opt_name) {
            if (v.empty())
                v.assign(cohort_num, default_val);
            else if (v.size() == 1)
                v.assign(cohort_num, v[0]);
            else if (v.size() != cohort_num)
                throw CLI::ValidationError(opt_name + " must provide 1 value or 1 value per cohort (" + to_string(cohort_num) + " cohorts)");
        };

        expand_or_validate(params.diff_freq_list, 0.2, "--diff-freq");
        expand_or_validate(params.maf_list, 0, "--maf");
        expand_or_validate(params.maf_sumstat_list, 0, "--maf-sumstat");
        expand_or_validate(params.missingness_list, 1, "--geno");
    });

    try { 
        app.parse(argc, argv); 
    } catch (const CLI::CallForHelp &e) {
        LOGGER.i(app.help());
        exit(0);
    } catch (const CLI::ParseError &e) {
        LOGGER.i(app.help());
        LOGGER.e(e.what());
    }

    params.if_joint_mode = joint_flag->count() > 0;
    params.if_cond_mode = cond_option->count() > 0;
    params.if_LD_mode = ld->count() > 0;

    params.window_size = (params.window_kb < 0 ? INT_MAX : params.window_kb * 1000);
    LOGGER.open(params.output_name + ".log");

    // check chromosome list across cohorts
    for (size_t i = 0; i < cohort_num; i++) {
        vector<int> cohort_chr_list, temp;

        // read all chromosomes in the .bim file
        string filename = params.bfile_list[i] + ".bim";
        ifstream sFile(filename.c_str());
        if (!sFile) LOGGER.e("Cannot open file [" + filename + "] to read");

        string str_buf;
        int chr_buf = -1;

        while (getline(sFile, str_buf)) {
            if (str_buf.empty()) continue;

            const char* pt = str_buf.c_str();
            parse_num(pt, chr_buf);
            cohort_chr_list.push_back(chr_buf);
        }

        sFile.close();

        sort(cohort_chr_list.begin(), cohort_chr_list.end());
        set_intersection(params.chr_list.begin(), params.chr_list.end(), cohort_chr_list.begin(), cohort_chr_list.end(), back_inserter(temp));
        params.chr_list.swap(temp);
    }

    if (params.chr_list.empty())
        LOGGER.e("No common chromosome found across all cohorts");

    if (params.curr_chr != -1) {
        if (find(params.chr_list.begin(), params.chr_list.end(), params.curr_chr) == params.chr_list.end()) {
            LOGGER.i("User specified chromosome " + to_string(params.curr_chr) + " with --chr option");
            LOGGER.e("No valid SNPs found on chromosome " + to_string(params.curr_chr));
        }
        
        params.chr_list.assign(1, params.curr_chr);
    }

    // set number of threads
    #if HAS_OPENMP 
        omp_set_num_threads(params.thread_num);
    #else
        if (params.thread_num > 1) {
            LOGGER.w("OpenMP is not available, running in single-thread mode");
            params.thread_num = 1;
        }
    #endif

    // output_user_hyperparameters
    LOGGER << "\n========== Manc-COJO CONFIGURATION ==========" << endl
            << "Thread Number: " << params.thread_num << endl;

    if (params.if_joint_mode)
        LOGGER << "Program Mode: Joint analysis only" << endl;
    else if (params.if_cond_mode)
        LOGGER << "Program Mode: Conditional analysis only" << endl;
    else
        LOGGER << "Program Mode: Stepwise iterative selection" << endl 
                << "Selection method: " + params.slct_mode << endl
                << "p-value threshold: " << params.p_value << endl;

    LOGGER << "Effect size estimation method: " << params.effect_size_mode << endl;
    
    if (params.if_LD_mode)
        LOGGER << "Input format: PLINK .ld files" << endl << endl;
    else
        LOGGER << "Input format: PLINK .bed files" << endl << endl;

    if (params.window_kb < 0)
        LOGGER << "No SNP position window applied" << endl;
    else
        LOGGER << "SNP position window (+/-): " << params.window_kb << " Kb" << endl;

    auto list_or_single = [](const vector<double>& v) -> string {
        auto fmt = [](double x) { ostringstream ss; ss << x; return ss.str(); };
        bool all_same = std::all_of(v.begin() + 1, v.end(), [&](double x){ return x == v[0]; });
        if (all_same) return fmt(v[0]);
        string s = "[";
        for (size_t i = 0; i < v.size(); i++) s += (i ? ", " : "") + fmt(v[i]);
        return s + "]";
    };

    LOGGER << "Collinearity threshold: " << params.collinear << endl
            << "Sumstat Minor Allele Frequency threshold: " << list_or_single(params.maf_sumstat_list) << endl
            << "Genotype Minor Allele Frequency threshold: " << list_or_single(params.maf_list) << endl
            << "Genotype missingness threshold: " << list_or_single(params.missingness_list) << endl
            << "Frequency difference threshold between sumstat and genotype: " << list_or_single(params.diff_freq_list) << endl;

    if (*R2_option)
        LOGGER << "R2 incremental threshold: " << params.R2_threshold << endl;
    if (*R2back_option)
        LOGGER << "R2 incremental threshold (backward): " << params.R2back_threshold << endl;
    if (params.if_var_from_ld)
        LOGGER << "Genotypic variance will be estimated from LD reference instead of 2pq" << endl;

    LOGGER << "=============================================" << endl << endl;

    auto load_snp_list = [&](const vector<string>& opts, vector<string>& snps,
                              const string& err_msg, const string& log_msg) {
        if (opts.empty() && snps.empty()) return;
        skim_SNP(opts, snps);
        if (snps.empty()) LOGGER.e(err_msg);
        LOGGER.i(log_msg, snps.size());
    };

    load_snp_list(params.extract_options, params.extract_SNPs,
        "Please provide valid SNPs for --extract and --extract-snp", "user-specified SNPs to include");
    load_snp_list(params.exclude_options, params.exclude_SNPs,
        "Please provide valid SNPs for --exclude and --exclude-snp", "user-specified SNPs to exclude");
    load_snp_list(params.fix_options, params.fix_SNPs,
        "Please provide valid SNPs for --fix and --fix-snp", "user-specified fixed SNPs");
    load_snp_list(params.cojo_cond_options, params.cojo_cond_SNPs,
        "Please provide valid SNPs for --cojo-cond", "user-specified conditional SNPs");
}
