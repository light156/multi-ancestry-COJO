#include "include/macojo.h"


// Parse CLI options, validate combinations, and populate HyperParams.
void set_read_process_output_options(int argc, char** argv)
{   
    HyperParams& params = get_params();
    CLI::App app;

    app.description(
        "\nMulti-ancestry conditional and joint analysis of GWAS summary statistics\n"
        "https://github.com/light156/multi-ancestry-COJO\n"
        "\nTIPS:\n"
        "  1. The usage is largely consistent with the original GCTA COJO, while extended to multiple cohorts and LD inputs. \n"
        "  2. Please make sure the paths are paired in --bfile/--ld and --cojo-file.\n"
        "     (e.g., --bfile xx1.sumstat xx2.sumstat xx3.sumstat --cojo-file yy1 yy2 yy3)\n"
        "  3. For using genotype data (--bfile), please provide (.bed .bim .fam) files.\n"
        "     For using LD matrix (--ld), please provide at least (.bim .ld) files.\n"
        "  4. The folders in the output file path (which is after --out) must exist.\n"
        "     (You can easily create them with 'mkdir -p <folder_path>')"
    );

    app.get_formatter()->column_width(75);
    app.get_formatter()->right_column_width(75);

    int cohort_num;

    // Main logic options
    string input_group = "Input Data Format Options";
    auto *bfile = app.add_option("--bfile", params.bfile_list, "PLINK binary file prefix for each cohort [.bim .bed .fam]")->group(input_group);
    auto *ld =  app.add_option("--ld", params.bfile_list, "PLINK LD file prefix for each cohort [.bim .ld (.frq)]")->group(input_group);
    
    string mode_group = "Algorithm Options/Flags";
    auto *slct_mode_option = app.add_option("--slct-mode", params.slct_mode, "Iterative SNP selection method")
        ->check(CLI::IsMember({"GCTA", "removeNA", "imputeNA"}))->default_val("GCTA")->group(mode_group);
    app.add_option("--effect-size-mode", params.effect_size_mode, "Effect size estimation method")
        ->check(CLI::IsMember({"GCTA", "removeNA", "imputeNA"}))->default_val("GCTA")->group(mode_group);

    // options same as GCTA-COJO
    string cojo_group = "Original GCTA COJO Options/Flags";
    app.add_option("--cojo-file", params.cojo_file_list, "GWAS summary statistics file for each cohort (GCTA-COJO format)")->required()->group(cojo_group);
    app.add_option("--out", params.output_name, "Output file path")->required()->group(cojo_group);
    auto *slct_flag = app.add_flag("--cojo-slct", "Stepwise iterative selection of independently associated SNPs")->group(cojo_group);
    auto *joint_flag = app.add_flag("--cojo-joint", "Only calculate joint effect for provided SNPs and exit")->group(cojo_group);
    auto *cond_option = app.add_option("--cojo-cond", params.cojo_cond_options, "Only calculate conditional effect for provided SNPs and exit")->group(cojo_group);
    app.add_option("--cojo-wind", params.window_kb, "SNP position window in Kb (-1 for no windows)")->default_val(10000)->group(cojo_group);
    auto *p_option = app.add_option("--cojo-p", params.p_value, "Significance threshold for SNP selection")->default_val(5e-8)->check(CLI::Range(0.0, 1.0))->group(cojo_group);
    app.add_option("--cojo-collinear", params.collinear, "Colinearity threshold")->default_val(0.9)->check(CLI::Range(0.0, 0.9999))->group(cojo_group);
    app.add_option("--diff-freq", params.diff_freq, "Frequency difference threshold between sumstat and genotype")->default_val(0.2)->check(CLI::Range(0.0, 1.0))->group(cojo_group);
    app.add_option("--maf", params.maf, "Minor Allele Frequency threshold")->default_val(0)->check(CLI::Range(0.0, 0.5))->group(cojo_group);
    
    string data_group = "Original GCTA Data Options";
    app.add_option("--geno", params.missingness, "Missingness threshold in PLINK .bed files (1 for no threshold)")->default_val(1.0)->check(CLI::Range(0.0, 1.0))->group("");
    app.add_option("--extract", params.extract_options, "File path of SNPs to be included")->group(data_group);
    app.add_option("--exclude", params.exclude_options, "File path of SNPs to be excluded")->group(data_group);
    app.add_option("--extract-snp", params.extract_SNPs, "A list of SNPs to be included")->group(data_group);
    app.add_option("--exclude-snp", params.exclude_SNPs, "A list of SNPs to be excluded")->group(data_group);
    auto *keep_option = app.add_option("--keep", params.keep_file_list, "File path of individuals to be included")->group(data_group);
    auto *remove_option = app.add_option("--remove", params.remove_file_list, "File path of individuals to be excluded")->group(data_group);
    app.add_option("--chr", params.curr_chr, "Only include SNPs on a specific chromosome")->group(data_group)->check(CLI::Range(1, 22));
    app.add_option("--thread-num", params.thread_num, "Number of threads to use")->default_val(1)->check(CLI::PositiveNumber)->group(data_group);
    app.add_option("--bed-block-mb", params.bed_block_mb, "BED read block size in MB")->default_val(64)->check(CLI::PositiveNumber)->group("");
    
    // multi-cohort specific options
    string manc_group = "Other Manc-COJO Options/Flags";
    auto *fix_option = app.add_option("--fix", params.fix_options, "File path of fixed SNPs for iterative selection")->group(manc_group);
    auto *fixSNP_option = app.add_option("--fix-snp", params.fix_SNPs, "A list of fixed SNPs for iterative selection")->group(manc_group);
    auto *R2_option = app.add_option("--R2", params.R2_threshold, "R2 incremental threshold (-1 for no threshold)")->default_val(-1)->group(manc_group);
    auto *R2back_option = app.add_option("--R2back", params.R2back_threshold, "R2 threshold for backward selection (-1 for no threshold)")->default_val(-1)->group(manc_group);
    auto *iter_option = app.add_option("--iter", params.max_iter_num, "Total iteration number")->default_val(10000)->check(CLI::PositiveNumber)->group(manc_group);
    app.add_flag("--freq-mode-and", params.if_freq_mode_and, "Use AND mode for frequency threshold across cohorts")->group(manc_group);
    app.add_flag("--adjust-infoscore", params.if_infoscore, "Adjust sumstat b")->group("");
    app.add_flag("--MDISA", params.if_MDISA, "Run single-ancestry analysis after multi-ancestry COJO")->group(manc_group);
    app.add_flag("--output-all", params.if_output_all, "Save all .cma.cojo .jma.cojo .ldr.cojo and .badsnps files")->group(manc_group);
    
    for (auto *opt : app.get_options())
        opt->multi_option_policy(CLI::MultiOptionPolicy::Throw);

    // Post-parse validation and derived parameter setup.
    app.callback([&]() {
        cohort_num = params.cojo_file_list.size();

        if ((slct_flag->count() > 0) + (joint_flag->count() > 0) + (cond_option->count() > 0) != 1)
            throw CLI::ValidationError("Exactly one of --cojo-slct, --cojo-joint, or --cojo-cond must be provided to specify the analysis mode");
        
        if ((bfile->count() > 0) == (ld->count() > 0))
            throw CLI::ValidationError("Either --bfile or --ld must be provided for all cohorts");

        if (*ld && (*keep_option || *remove_option))
            throw CLI::ValidationError("--keep/--remove cannot be used with --ld");

        if (params.bfile_list.size() != cohort_num)
            throw CLI::ValidationError("Filepath numbers after --bfile/--ld and --cojo-file must be the same");

        if (*keep_option && params.keep_file_list.size() != cohort_num)
            throw CLI::ValidationError("Filepath numbers after --cojo-file and --keep must be the same");

        if (*remove_option && params.remove_file_list.size() != cohort_num)
            throw CLI::ValidationError("Filepath numbers after --cojo-file and --remove must be the same");

        if ((*joint_flag || *cond_option) && (*slct_mode_option || *p_option || *fix_option || *fixSNP_option || *iter_option || *R2_option || *R2back_option || params.if_MDISA))
            throw CLI::ValidationError("--cojo-joint and --cojo-cond cannot be used with selection-related options such as --slct-mode/--cojo-p/--fix/--fix-snp/--iter/--R2/--R2back/--MDISA, because it makes no sense");
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
    params.iter_collinear_threshold = 1.0 / (1.0 - params.collinear);
    LOGGER.open(params.output_name + ".log");

    if (params.keep_file_list.empty())
        params.keep_file_list.assign(cohort_num, "");

    if (params.remove_file_list.empty())
        params.remove_file_list.assign(cohort_num, "");

    // check chromosome list across cohorts
    for (int i = 0; i < cohort_num; i++) {
        vector<int> cohort_chr_list, temp;

        // read all chromosomes in the .bim file
        string filename = params.bfile_list[i] + ".bim";
        ifstream sFile(filename.c_str());
        if (!sFile) LOGGER.e("Cannot open file [" + filename + "] to read");

        string str_buf;
        int chr_buf;

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

    if (params.chr_list.size() == 0)
        LOGGER.e("No common chromosome found across all cohorts");

    if (params.curr_chr != -1) {
        if (find(params.chr_list.begin(), params.chr_list.end(), params.curr_chr) == params.chr_list.end()) {
            LOGGER.i("User specified chromosome " + to_string(params.curr_chr) + " with --chr option");
            LOGGER.e("No valid SNPs found on chromosome " + to_string(params.curr_chr));
        }
        
        params.chr_list = {params.curr_chr};
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
    LOGGER << "\n=========== MACOJO CONFIGURATION ===========" << endl
            << "Thread Number: " << params.thread_num << endl;

    if (params.if_joint_mode)
        LOGGER << "Program Mode: Joint analysis only" << endl;
    else if (params.if_cond_mode)
        LOGGER << "Program Mode: Conditional analysis only" << endl;
    else
        LOGGER << "Program Mode: Stepwise iterative selection" << endl 
                << "Selection method: " + params.slct_mode << endl
                << "Effect size estimation method: " << params.effect_size_mode << endl;
    
    if (params.if_LD_mode)
        LOGGER << "Input format: PLINK .ld files" << endl << endl;
    else
        LOGGER << "Input format: PLINK .bed files" << endl << endl;

    if (params.window_kb < 0)
        LOGGER << "No SNP position window applied" << endl;
    else
        LOGGER << "SNP position window (+/-): " << params.window_kb << " Kb" << endl;

    LOGGER << "p-value Threshold: " << params.p_value << endl
            << "Collinearity threshold: " << params.collinear << endl
            << "Sumstat Minor Allele Frequency threshold: " << params.maf << endl
            // << "Genotype Missingness threshold: " << params.missingness << endl
            << "Frequency difference threshold between sumstat and genotype: " << params.diff_freq << endl
            << (*R2_option ? "R2 incremental threshold: " + to_string(params.R2_threshold) + "\n" : "")
            << (*R2back_option ? "R2 incremental threshold (backward): " + to_string(params.R2back_threshold) + "\n" : "")
            << (cohort_num > 1 ? (params.if_freq_mode_and ? "SNP frequency mode: AND\n" : "SNP frequency mode: OR\n") : "")
            << (cohort_num > 1 && params.if_MDISA ? "Run MDISA after MACOJO\n" : "")
            << "===========================================" << endl << endl;

    if (!params.extract_options.empty() || !params.extract_SNPs.empty()) {
        skim_SNP(params.extract_options, params.extract_SNPs);
        if (params.extract_SNPs.size() == 0)
            LOGGER.e("Please provide valid SNPs for --extract and --extract-snp");

        LOGGER.i("user-specified SNPs to include", params.extract_SNPs.size());
    }

    if (!params.exclude_options.empty() || !params.exclude_SNPs.empty()) {
        skim_SNP(params.exclude_options, params.exclude_SNPs);
        if (params.exclude_SNPs.size() == 0)
            LOGGER.e("Please provide valid SNPs for --exclude and --exclude-snp");
            
        LOGGER.i("user-specified SNPs to exclude", params.exclude_SNPs.size());
    }

    if (!params.fix_options.empty() || !params.fix_SNPs.empty()) {
        skim_SNP(params.fix_options, params.fix_SNPs);
        if (params.fix_SNPs.size() == 0)
            LOGGER.e("Please provide valid SNPs for --fix and --fix-snp");

        LOGGER.i("user-specified fixed SNPs", params.fix_SNPs.size());
    }

    if (!params.cojo_cond_options.empty()) {
        skim_SNP(params.cojo_cond_options, params.cojo_cond_SNPs);
        if (params.cojo_cond_SNPs.size() == 0)
            LOGGER.e("Please provide valid SNPs for --cojo-cond");

        LOGGER.i("user-specified conditional SNPs", params.cojo_cond_SNPs.size());
    }
}
