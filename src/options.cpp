#include "macojo.h"


int MACOJO::set_read_process_output_options(int argc, char** argv)
{   
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

    int thread_num, cohort_num;

    // Main logic options
    string input_group = "Input Data Format Options";
    auto *bfile = app.add_option("--bfile", params.bfile_list, "PLINK binary file prefix for each cohort [.bim .bed .fam]")->group(input_group);
    auto *ld =  app.add_option("--ld", params.bfile_list, "PLINK LD file prefix for each cohort [.bim .ld (.frq)]")->group(input_group);
    // app.add_option("--bgen", params.if_bgen_mode, "PLINK BGEN file prefix for each cohort (.bim .bgen .fam)")->group(input_group);
    
    string mode_group = "Algorithm Options/Flags";
    auto *slct_mode_option = app.add_option("--slct-mode", params.slct_mode, "Iterative SNP selection method")
        ->check(CLI::IsMember({"GCTA", "removeNA", "imputeNA"}))->default_val("GCTA")->group(mode_group);
    app.add_option("--effect-size-mode", params.effect_size_mode, "Effect size estimation method")
        ->check(CLI::IsMember({"GCTA", "removeNA", "imputeNA"}))->default_val("GCTA")->group(mode_group);

    // options same as GCTA-COJO
    string cojo_group = "Original GCTA COJO Options/Flags";
    auto *slct_flag = app.add_flag("--cojo-slct", "Stepwise iterative selection of independently associated SNPs")->group(cojo_group);
    auto *joint_flag = app.add_flag("--cojo-joint", "Only calculate joint effect for provided SNPs and exit")->group(cojo_group);
    auto *cond_option = app.add_option("--cojo-cond", params.cond_file, "Only calculate conditional effect for provided SNPs and exit")->check(CLI::ExistingFile)->group(cojo_group);
    app.add_option("--cojo-file", params.cojo_file_list, "GWAS summary-level statistics file for each cohort")->required()->check(CLI::ExistingFile)->group(cojo_group);
    app.add_option("--out", params.output_name, "Output file path")->required()->group(cojo_group);
    app.add_option("--cojo-wind", params.window_kb, "SNP position window in Kb (-1 for no windows)")->default_val(10000)->group(cojo_group);
    auto *p_option = app.add_option("--cojo-p", params.p, "Significance threshold for SNP selection")->default_val(5e-8)->check(CLI::Range(0.0, 1.0))->group(cojo_group);
    app.add_option("--cojo-collinear", params.collinear, "Colinearity threshold (0-0.9999)")->default_val(0.9)->check(CLI::Range(0.0, 0.9999))->group(cojo_group);
    app.add_option("--diff-freq", params.diff_freq, "Frequency difference threshold between sumstat and PLINK (0-1)")->default_val(0.2)->check(CLI::PositiveNumber)->group(cojo_group);

    string data_group = "Original GCTA Data Options";
    app.add_option("--maf", params.maf, "Minor Allele Frequency threshold")->default_val(0.01)->check(CLI::Range(1e-5, 0.5))->group(data_group);
    app.add_option("--geno", params.missingness, "Missingness threshold in PLINK .bed files(-1 for no threshold)")->default_val(1.0)->check(CLI::Range(0.0, 1.0))->group(data_group);
    auto *extract_option = app.add_option("--extract", params.extract_file, "File path of SNPs to be included")->check(CLI::ExistingFile)->group(data_group);
    auto *exclude_option = app.add_option("--exclude", params.exclude_file, "File path of SNPs to be excluded")->check(CLI::ExistingFile)->group(data_group);
    auto *keep_option = app.add_option("--keep", params.keep_file_list, "File path of individuals to be included")->check(CLI::ExistingFile)->group(data_group);
    auto *remove_option = app.add_option("--remove", params.remove_file_list, "File path of individuals to be excluded")->check(CLI::ExistingFile)->group(data_group);
    
    joint_flag->needs(extract_option);
    
    // multi-cohort specific options
    string manc_group = "Other Manc-COJO Options/Flags";
    auto *fixedSNP_option = app.add_option("--fixed", params.fixedSNP_file, "File path of fixed candidate SNPs")->check(CLI::ExistingFile)->group(manc_group);
    auto *R2_option = app.add_option("--R2", params.R2_threshold, "R2 incremental threshold (-1 for no threshold)")->default_val(-1)->group(manc_group);
    auto *R2back_option = app.add_option("--R2back", params.R2back_threshold, "R2 threshold for backward selection (-1 for no threshold)")->default_val(-1)->group(manc_group);
    auto *iter_option = app.add_option("--iter", params.max_iter_num, "Total iteration number")->default_val(10000)->check(CLI::PositiveNumber)->group(manc_group);
    app.add_option("--thread-num", thread_num, "Number of threads to use")->default_val(1)->check(CLI::Range(1,10))->group(manc_group);
    app.add_flag("--freq-mode-and", params.if_freq_mode_and, "Use AND mode for frequency threshold across cohorts")->group(manc_group);
    app.add_flag("--MDISA", params.if_MDISA, "Run single-ancestry analysis after multi-ancestry COJO")->group(manc_group);
    app.add_flag("--output-all", params.if_output_all, "Save all .cma.cojo, .jma.cojo and .ldr.cojo results to file")->group(manc_group);

    for (auto *opt : app.get_options())
        opt->multi_option_policy(CLI::MultiOptionPolicy::Throw);

    app.callback([&]() {
        cohort_num = params.cojo_file_list.size();

        if (slct_flag->count() + joint_flag->count() + cond_option->count() != 1)
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

        if (params.if_MDISA && cohort_num == 1)
            throw CLI::ValidationError("--MDISA are invalid options with single cohort");

        if ((*joint_flag || *cond_option) && (*slct_mode_option || *p_option || *fixedSNP_option || *iter_option || *R2_option || *R2back_option || params.if_MDISA))
            throw CLI::ValidationError("--cojo-joint and --cojo-cond cannot be used with --slct-mode/--cojo-p/--fixed/--iter/--R2/--R2back/--MDISA, because it makes no sense");

        if (params.slct_mode == "gcta" && (*R2_option || *R2back_option))
            throw CLI::ValidationError("--slct-mode GCTA cannot be used with --R2/--R2back, because GCTA-COJO does not check R2 increment");
    });

    try { 
        app.parse(argc, argv); 
    } catch (const CLI::CallForHelp &e) {
        LOGGER << app.help() << endl;
        return -1;
    } catch (const CLI::ParseError &e) {
        LOGGER << app.help() << endl;
        LOGGER.e(0, e.what());
    }

    params.if_joint_mode = joint_flag->count() > 0;
    params.if_cond_mode = cond_option->count() > 0;
    params.if_LD_mode = ld->count() > 0;

    params.window_size = (params.window_kb < 0 ? INT_MAX : params.window_kb * 1000);
    params.iter_collinear_threshold = 1.0 / (1.0 - params.collinear);
    LOGGER.open(params.output_name + ".log");

    for (size_t i = 0; i < cohort_num; i++)
        cohorts.emplace_back(params, shared, i);

    // output_user_hyperparameters
    LOGGER << "\n=========== MACOJO CONFIGURATION ===========" << endl
            << (params.if_joint_mode ? "Program Mode: Joint analysis only\n" : 
                (params.if_cond_mode ? "Program Mode: Conditional analysis only\n" : 
                    "Program Mode: Stepwise iterative selection\nSelection method: " + params.slct_mode + "\n"))
            << "Effect size estimation method: " << params.effect_size_mode << "\n"
            << (params.if_LD_mode ? "Input format: PLINK .ld files\n" : "Input format: PLINK .bed files\n")
            << "\n"
            << "p-value Threshold: " << params.p << endl
            << "Collinearity threshold: " << params.collinear << endl
            << "SNP position window (+/-): " << params.window_kb << "Kb" << endl
            << "Sumstat Minor Allele Frequency threshold: " << params.maf << endl
            << "Genotype Missingness threshold: " << params.missingness << endl
            << "Frequency difference threshold between sumstat and PLINK: " << params.diff_freq << endl
            << (*R2_option ? "R2 incremental threshold: " + to_string(params.R2_threshold) + "\n" : "")
            << (*R2back_option ? "R2 incremental threshold (backward): " + to_string(params.R2back_threshold) + "\n" : "")
            << (cohorts.size() > 1 ? 
                (params.if_freq_mode_and ? "SNP frequency mode: AND\n" : "SNP frequency mode: OR\n") : "")
            << (params.if_MDISA ? "Run MDISA after MACOJO\n" : "")
            << "===========================================" << endl << endl;

    #if HAS_OPENMP
        omp_set_num_threads(thread_num);
    #else
        if (thread_num > 1)
            LOGGER.w(0, "OpenMP is not available, running in single-threaded mode");
    #endif

    return 0;
}
