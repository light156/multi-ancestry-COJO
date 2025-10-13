#include "macojo.h"
#include "CLI11.hpp"


int main(int argc, char** argv) 
{     
    MACOJO macojo;
    string extract_file, fixedSNP_file, savename;
    vector<string> holder;

    CLI::App app;
    app.get_formatter()->column_width(60);
    app.description(
        "Multi-ancestry conditional and joint analysis of GWAS summary statistics\n"
        "TIPS:\n"
        "  1. The first parameter indicates the number of cohorts, which should be a positive integer.\n"
        "  2. The next 2*N parameters are file paths of GWAS summary statistics and PLINK files for N cohorts, respectively.\n"
        "  3. All options must be placed after the above parameters.\n"
        "  4. The folders in the output file path (which is after -o or --out) must exist.\n"
        "     (You can easily create them with 'mkdir -p <folder_path>')"
    );

    app.add_option("holder", holder, "Cohort number followed by sumstat and PLINK file paths")->take_all()->group("");
    app.add_option("-o, --out", savename, "Output file path")->required();
    app.add_option("-e, --extract", extract_file, "File path of user-given SNPs")->check(CLI::ExistingFile);
    app.add_option("-f, --fixedSNP", fixedSNP_file, "File path of fixed candidate SNPs")->check(CLI::ExistingFile);
    app.add_option("-w, --window", macojo.window_mb, "SNP position window in Mb (-1 for no windows)")->default_val(10.0);
    app.add_option("-c, --colinear", macojo.colinear_threshold, "Colinearity threshold (0-1)")->default_val(0.9)->check(CLI::Range(0.0, 1.0));
    app.add_option("--R2", macojo.R2_incremental_threshold, "R2 incremental threshold (-1 for no threshold)")->default_val(-1);
    app.add_option("--R2back", macojo.R2_incremental_threshold_backwards, "R2 threshold for backward selection (-1 for no threshold)")->default_val(-1);
    app.add_option("--freq", macojo.freq_threshold, "Frequency threshold (0-0.5)")->default_val(0.01)->check(CLI::Range(0.0, 0.5));
    app.add_option("--iter", macojo.max_iter_num, "Total iteration number")->default_val(10000)->check(CLI::PositiveNumber);

    for (auto *opt : app.get_options())
        opt->multi_option_policy(CLI::MultiOptionPolicy::Throw);

    app.add_flag("--freq_mode_and", macojo.if_freq_mode_and, "Use AND mode for frequency threshold");
    app.add_flag("--LD", macojo.if_LD_mode, "Use .ld files for calculation");
    app.add_flag("--cojo-joint", macojo.if_cojo_joint, "Only output for provided fixed candidate SNPs and exit");
    app.add_flag("--skip_MDISA", macojo.if_skip_MDISA, "Disable MDISA after COJO");
    app.add_flag("--keep_NA", macojo.if_keep_NA, "Do not fill NA with mean values");
    app.add_flag("--fast_inv", macojo.if_fast_inv, "Incremental matrix inverse based on last iteration");
    app.add_flag("--gcta", macojo.if_gcta_COJO, "Use GCTA-COJO model selection criteria");

    if (argc < 2 || atoi(argv[1]) <= 0) {
        LOGGER << app.help() << endl;
        return -1;
    }  

    int cohort_num = atoi(argv[1]);
    macojo.cohorts.resize(cohort_num);
    if (argc < cohort_num * 2 + 4) {
        LOGGER << "Too few parameters, please check how to use this program." << endl;
        LOGGER << app.help() << endl;
        return -1;
    }  

    CLI11_PARSE(app, argc, argv);
    LOGGER.open(savename + ".log");
    LOGGER << setprecision(12);

    // double tStart = clock();
    double tStart = omp_get_wtime();
    
    if (!extract_file.empty()) {
        macojo.read_SNP_only(extract_file, macojo.all_SNP);
        LOGGER.i(0, "The user has initialized SNPs for analysis");
    }

    if (!fixedSNP_file.empty()) {
        macojo.read_SNP_only(fixedSNP_file, macojo.fixed_candidate_SNP);
        LOGGER << "The user has provided fixed candidate SNPs, which will not be removed during calculation" << endl;
    }

    macojo.output_user_hyperparameters();
    macojo.read_cojo_PLINK_files(argv+2, cohort_num);
    macojo.entry_function(savename);

    // LOGGER << "Total running time: " << fixed << setprecision(2) << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "Total running time: " << fixed << setprecision(2) << omp_get_wtime() - tStart << " seconds" << endl;
    LOGGER.close();

    return 0;
}