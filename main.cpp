#include "macojo.h"
#include "CLI11.hpp"


int main(int argc, char** argv) 
{    
    LOGGER.ts("main");

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
    app.add_option("-w, --window", macojo.params.window_mb, "SNP position window in Mb (-1 for no windows)")->default_val(10.0);
    
    app.add_option("--cojo-p", macojo.params.threshold, "Significance threshold for SNP selection")->default_val(5e-8)->check(CLI::Range(0.0, 1.0));
    app.add_option("--cojo-collinear", macojo.params.colinear_threshold, "Colinearity threshold (0-0.999)")->default_val(0.9)->check(CLI::Range(0.0, 0.999));
    app.add_option("--R2", macojo.params.R2_incremental_threshold, "R2 incremental threshold (-1 for no threshold)")->default_val(-1);
    app.add_option("--R2back", macojo.params.R2_incremental_threshold_backwards, "R2 threshold for backward selection (-1 for no threshold)")->default_val(-1);
    app.add_option("--freq", macojo.params.freq_threshold, "Frequency threshold (0-0.5)")->default_val(0.01)->check(CLI::Range(0.0, 0.5));
    app.add_option("--diff-freq", macojo.params.freq_diff_threshold, "Frequency difference threshold between sumstat and PLINK (0-1)")->default_val(0.2)->check(CLI::Range(0.0, 1.0));
    app.add_option("--iter", macojo.params.max_iter_num, "Total iteration number")->default_val(10000)->check(CLI::PositiveNumber);

    for (auto *opt : app.get_options())
        opt->multi_option_policy(CLI::MultiOptionPolicy::Throw);

    app.add_flag("--freq-mode-and", macojo.params.if_freq_mode_and, "Use AND mode for frequency threshold");
    app.add_flag("--LD", macojo.params.if_LD_mode, "Use .ld files for calculation");
    app.add_flag("--cojo-joint", macojo.params.if_cojo_joint, "Only output for provided fixed candidate SNPs and exit");
    app.add_flag("--MDISA", macojo.params.if_MDISA, "Run MDISA after MACOJO");
    app.add_flag("--keep-NA", macojo.params.if_keep_NA, "Do not fill NA with mean values");
    app.add_flag("--gcta", macojo.params.if_gcta_COJO, "Use GCTA-COJO model selection criteria");
    // app.add_flag("--fast-inv", macojo.params.if_fast_inv, "Incremental matrix inverse based on last iteration");

    if (argc < 2 || atoi(argv[1]) <= 0) {
        LOGGER << app.help() << endl;
        return -1;
    }  

    int cohort_num = atoi(argv[1]);
    if (argc < cohort_num * 2 + 4) {
        LOGGER << "Too few parameters, please check how to use this program." << endl;
        LOGGER << app.help() << endl;
        return -1;
    }  

    CLI11_PARSE(app, argc, argv);
    macojo.params.window_size = (macojo.params.window_mb < 0 ? INT_MAX : macojo.params.window_mb * 1e6);
    macojo.params.iter_colinear_threshold = 1.0 / (1.0 - macojo.params.colinear_threshold);
    LOGGER.open(savename + ".log");
    LOGGER << setprecision(12);

    for (int i = 0; i < cohort_num; i++)
        macojo.cohorts.emplace_back(macojo.params, macojo.shared);

    if (!fixedSNP_file.empty()) {    
        macojo.read_SNP_only(fixedSNP_file, macojo.fixed_candidate_SNP);
        LOGGER.i(0, "fixed candidate SNPs provided for analysis", to_string(macojo.fixed_candidate_SNP.size()));
    }

    macojo.output_user_hyperparameters();
    macojo.read_cojo_PLINK_files(argv+2, cohort_num, extract_file);
    macojo.entry_function(savename);

    #ifdef __linux__
        float vmem = roundf(1000.0 * getVMPeakKB() / 1024/1024) / 1000.0; 
        float mem = roundf(1000.0 * getMemPeakKB() / 1024/1024) / 1000.0;     
        LOGGER << setprecision(3) << "Peak memory: " << mem << " GB; Virtual memory: " << vmem << " GB." << endl;
    #endif

    float duration = LOGGER.tp("main");
    int hours = (int) duration / 3600;
    string time_str = (hours == 0) ? "" : (to_string(hours) + " hour" + ((hours == 1) ? " ": "s "));
    int mins = (int) (duration - 3600 * hours) / 60;
    time_str += (mins == 0) ? "" : (to_string(mins) + " minute" + ((mins == 1) ? " ": "s "));
    float seconds = duration - 3600 * hours - 60 * mins;
    time_str = time_str + to_string(seconds) + " seconds";
    LOGGER << "Overall computational time: " << time_str << endl;
    
    LOGGER.close();
    return 0;
}