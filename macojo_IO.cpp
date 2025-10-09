#include "macojo.h"


void MACOJO::output_results_to_file(string filepath) 
{   
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    jmaCOJO << "SNP" << "\t" 
            << "A1" << "\t" 
            << "A2" << "\t";

    for (int n : current_calculation_list)
        jmaCOJO << "freq." << to_string(n+1) << "\t" 
                << "b." << to_string(n+1) << "\t" 
                << "se." << to_string(n+1) << "\t" 
                << "p." << to_string(n+1) << "\t" 
                << "N." << to_string(n+1) << "\t"
                << "bJ." << to_string(n+1) << "\t" 
                << "seJ." << to_string(n+1) << "\t" 
                << "zJ." << to_string(n+1) << "\t" 
                << "pJ." << to_string(n+1) << "\t";

    ArrayXd bJma, seJma, zJma, pJma;

    if (current_calculation_list.size() > 1) {
        jmaCOJO << "bJ.ma" << "\t" 
                << "seJ.ma" << "\t" 
                << "zJ.ma" << "\t" 
                << "pJ.ma" << "\t";
                
        // calculate output sumstats    
        sumstat_new_model_joint.setZero(candidate_SNP.size(), 4);

        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            sumstat_new_model_joint.col(0) += c.output_b / c.output_se2;
            sumstat_new_model_joint.col(1) += 1 / c.output_se2;
        }

        sumstat_new_model_joint.col(1) = 1 / sumstat_new_model_joint.col(1);
        sumstat_new_model_joint.col(0) = sumstat_new_model_joint.col(0) * sumstat_new_model_joint.col(1);

        bJma = sumstat_new_model_joint.col(0);
        seJma = sqrt(sumstat_new_model_joint.col(1));
        zJma = bJma / seJma;
        pJma = erfc(abs(zJma)/sqrt(2));
    } 
    
    jmaCOJO << "\n";

    int index = 0;
    map<string, int> SNP_order_pair;
 
    // Inserting element in pair vector to keep track of previous indexes
    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, index++) 
        SNP_order_pair.insert(make_pair(final_commonSNP[*iter], index));
    
    jmaCOJO.precision(12);

    for (auto iter = SNP_order_pair.begin(); iter != SNP_order_pair.end(); iter++) {
        index = iter->second;
        jmaCOJO << iter->first << "\t" 
                << A1_ref[final_commonSNP_index[candidate_SNP[index]]] << "\t" 
                << A2_ref[final_commonSNP_index[candidate_SNP[index]]] << "\t";

        for (int n : current_calculation_list) {
            auto &c = cohorts[n];
            
            // cols 0:b, 1:se2, 2:p, 3:freq, 4:N
            jmaCOJO << c.sumstat_candidate(index, 3) << "\t" 
                    << c.sumstat_candidate(index, 0) << "\t" 
                    << sqrt(c.sumstat_candidate(index, 1)) << "\t" 
                    << c.sumstat_candidate(index, 2) << "\t" 
                    << c.sumstat_candidate(index, 4) << "\t"
                    << c.output_b(index) << "\t"
                    << sqrt(c.output_se2(index)) << "\t"
                    << c.output_b(index) / sqrt(c.output_se2(index)) << "\t"
                    << erfc(abs(c.output_b(index)) / sqrt(c.output_se2(index)) / sqrt(2)) << "\t";
        }

        if (current_calculation_list.size() > 1) {
            jmaCOJO << bJma(index) << "\t" 
                    << seJma(index) << "\t" 
                    << zJma(index) << "\t" 
                    << pJma(index) << "\t";
        } 
        
        jmaCOJO << "\n";
    }

    jmaCOJO.clear();
    jmaCOJO.close();

    map<string, int>().swap(SNP_order_pair);
    LOGGER.i(0, "Results saved into [" + filepath + "]\n", "Finished");
}


void MACOJO::read_user_hyperparameters(int argc, char** argv) 
{
    int temp_num = 0;

    while (temp_num < argc) { 
        if (strcmp(argv[temp_num], "-iter_num") == 0 && temp_num+1 < argc) {
            max_iter_num = atoi(argv[temp_num+1]);
            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "-extract") == 0 && temp_num+1 < argc) {
            read_SNP_only(argv[temp_num+1], all_SNP, false);
            LOGGER << "The user has initialized SNPs for analysis" << endl << endl;
            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "-fixedSNP") == 0 && temp_num+1 < argc) {
            read_SNP_only(argv[temp_num+1], fixed_candidate_SNP, false);
            LOGGER << "The user has provided fixed candidate SNPs, which will not be removed during calculation" << endl << endl;
            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "-colinear") == 0 && temp_num+1 < argc) {
            colinear_threshold = atof(argv[temp_num+1]);
            if (colinear_threshold < 0 || colinear_threshold >= 1) 
                LOGGER.e(0, "colinear_threshold should be between 0 and 1");

            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "-R2") == 0 && temp_num+1 < argc) {
            R2_incremental_threshold = atof(argv[temp_num+1]);
            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "-R2back") == 0 && temp_num+1 < argc) {
            R2_incremental_threshold_backwards = atof(argv[temp_num+1]);
            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "-window") == 0 && temp_num+1 < argc) {
            if (atof(argv[temp_num+1]) > 0)
                window_size = atof(argv[temp_num+1])*1e6;
            else
                window_size = INT_MAX;

            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "-freq") == 0 && temp_num+1 < argc) {    
            freq_threshold = atof(argv[temp_num+1]);
            if (freq_threshold < 0 || freq_threshold >= 0.5) 
                LOGGER.e(0, "freq threshold should be between 0 and 0.5");

            temp_num += 2;
        } 
        else if (strcmp(argv[temp_num], "--freq_mode_and") == 0) {
            if_freq_mode_or = false;
            temp_num += 1;
        } 
        else if (strcmp(argv[temp_num], "--LD") == 0) {
            if_LD_mode = true;
            temp_num += 1;
        }
        else if (strcmp(argv[temp_num], "--no_MDISA") == 0) {
            if_MDISA = false;
            temp_num += 1;
        } 
        else if (strcmp(argv[temp_num], "--no_fast_inv") == 0) {
            if_fast_inv = false;
            temp_num += 1;
        } 
        else if (strcmp(argv[temp_num], "--cojo-joint") == 0) {
            if_cojo_joint = true;
            temp_num += 1;
        } 
        else if (strcmp(argv[temp_num], "--fill_NA") == 0) {
            if_fill_NA = true;
            temp_num += 1;
        }
        else if (strcmp(argv[temp_num], "--gcta") == 0) {
            if_gcta_COJO = true;
            temp_num += 1;
        }
        else {
            show_tips();
            LOGGER.e(0, "Unknown option: " + string(argv[temp_num]));
        }
    }
    
    LOGGER << "Threshold: 5e-8" << endl;
    LOGGER << "Colinearity threshold: " << colinear_threshold << endl;
    LOGGER << "R2 incremental threshold: " << R2_incremental_threshold << endl;
    LOGGER << "R2 incremental threshold backwards: " << R2_incremental_threshold_backwards << endl;
    LOGGER << "SNP position window (+/-): " << window_size << endl;
    LOGGER << "SNP frequency threshold: " << freq_threshold << endl;
    
    if (if_freq_mode_or)
        LOGGER << "SNP frequency mode: OR" << endl;
    else
        LOGGER << "SNP frequency mode: AND" << endl;

    if (if_LD_mode)
        LOGGER << endl << "Use .ld files for calculation" << endl;

    if (!if_MDISA)
        LOGGER << endl << "Do not run MDISA after MACOJO" << endl;

    if (if_fill_NA)
        LOGGER << endl << "Fill NA with mean values for genotypes" << endl;

    for (auto &c : cohorts) {
        c.window_size = window_size;
        c.if_fast_inv = if_fast_inv;
        c.if_LD_mode = if_LD_mode;
        c.if_fill_NA = if_fill_NA;
        c.if_gcta_COJO = if_gcta_COJO;
        c.iter_colinear_threshold = 1 / (1 - colinear_threshold);
    }

    LOGGER << "--------------------------------" << endl << endl;
}


void MACOJO::show_tips() 
{
    cout << endl << "Multi-ancestry conditional and joint analysis of GWAS summary statistics" << endl;
    cout << "Usage: <program_path> <cohort_num> <sumstat_file1> <PLINK_file1> [<sumstat_file2> <PLINK_file2> ...] <output_file> [options]" << endl;
    
    cout << endl << "Examples:" << endl;
    cout << "Two cohorts: ./manc_cojo 2 sumstat1_path PLINK1_path sumstat2_path PLINK2_path output_path -colinear 0.9" << endl;
    cout << "One cohort: ./manc_cojo 1 sumstat_path PLINK_path output_path -fixedSNP user_file -window 10" << endl;

    cout << endl << "TIPS:" << endl;
    cout << "1. The first parameter indicates the number of cohorts, which should be a positive integer" << endl;
    cout << "2. Any folders in output_path will not be created, please make sure they exist" << endl;
    cout << "3. Please put all options after output_path" << endl << endl;

    cout << "Options:" << endl;   
    cout << "-extract: file path of user-given SNPs included for analysis" << endl;
    cout << "-fixedSNP: file path of fixed candidate SNPs, not removable during iterations" << endl;
    cout << "-colinear: colinear_threshold (default: 0.9)" << endl;
    cout << "-R2: R2_incremental_threshold (default: -1, no threshold)" << endl;
    cout << "-R2back: R2_incremental_threshold for backward selection (default: -1, no threshold)" << endl;
    cout << "-window (Mb): SNP position window +/- (default: 10 [which means +/-10Mb], set to -1 for no window" << endl;
    cout << "-freq: frequency threshold to exclude rare SNPs (default: 0.01)" << endl;
    cout << "--freq_mode_and: only keep SNPs that reach frequency threshold in sumstat files of all cohorts (default: at least one cohort)" << endl;   
    cout << "--LD: read .ld files instead of .bed files, please refer to GitHub for details" << endl;
    cout << "--no_MDISA: do not run MDISA after multi-ancestry COJO (can be ignored for single ancestry)" << endl;
    cout << "--fill_NA: fill NA with mean values for genotypes during inner product calculation (default: no filling)" << endl;
    cout << "--cojo-joint: only output for provided fixed candidate SNPs and exit" << endl;
    cout << "-iter_num: total iteration number (default: 10000)" << endl;
    cout << "--gcta: use GCTA-COJO model selection criteria (default: do not use)" << endl;
    
    // below are some trivial features, hidden in the usage tips 
    // cout << "--no_fast_inv: use normal matrix inverse (default: use fast inverse)" << endl;

    cout << endl;
}
