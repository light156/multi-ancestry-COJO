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
        if (strcmp(argv[temp_num], "-extract") == 0 && temp_num+1 < argc) {
            read_SNP_only(argv[temp_num+1], all_SNP, false);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-fixedSNP") == 0 && temp_num+1 < argc) {
            read_SNP_only(argv[temp_num+1], fixed_candidate_SNP, false);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-colinear") == 0 && temp_num+1 < argc) {
            colinear_threshold = atof(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-R2") == 0 && temp_num+1 < argc) {
            R2_incremental_threshold = atof(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-R2back") == 0 && temp_num+1 < argc) {
            R2_incremental_threshold_backwards = atof(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-iter_num") == 0 && temp_num+1 < argc) {
            max_iter_num = atoi(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-window") == 0 && temp_num+1 < argc && atoi(argv[temp_num+1]) == -1) {
            window_size = INT_MAX;
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-window") == 0 && temp_num+1 < argc && atof(argv[temp_num+1]) > 0) {
            window_size = atof(argv[temp_num+1])*1e6;
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "--no_fast_inv") == 0) {
            if_fast_inv = false;
            temp_num += 1;
        } else if (strcmp(argv[temp_num], "--no_MDISA") == 0) {
            if_MDISA = false;
            temp_num += 1;
        } else if (strcmp(argv[temp_num], "--cojo-joint") == 0) {
            if_cojo_joint = true;
            temp_num += 1;
        } else {
            LOGGER.e(0, "Unknown option: " + string(argv[temp_num]));
            show_tips_and_exit();
        }
    }

    colinear_threshold_sqrt = sqrt(colinear_threshold);
    iter_colinear_threshold = 1 / (1 - colinear_threshold);

    LOGGER << "Threshold: 5e-8" << endl;
    LOGGER << "Colinearity threshold: " << colinear_threshold << endl;
    LOGGER << "R2 incremental threshold: " << R2_incremental_threshold << endl;
    LOGGER << "R2 incremental threshold backwards: " << R2_incremental_threshold_backwards << endl;
    LOGGER << "Maximal iteration number: " << max_iter_num << endl;
    LOGGER << "SNP position window (+/-): " << window_size << endl;
    if (!if_fast_inv)
        LOGGER << "Use normal matrix inversion" << endl;
    
    if (!if_MDISA)
        LOGGER << "Do not run MDISA after MACOJO" << endl;

    LOGGER << "--------------------------------" << endl;
}


void MACOJO::show_tips_and_exit() 
{
    cout << endl << "Program usage: Joint meta analysis for two cohorts (--d) or a single cohort (--s)" << endl;
    cout << "Usage 1: program_path cohort_number cojoFile1_path PLINK1_path cojoFile2_path PLINK2_path result_save_path (+other options)" << endl;
    cout << "Usage 2: program_path cohort_number cojoFile_path PLINK_path result_save_path (+other options)" << endl << endl;

    cout << "TIPS:" << endl;
    cout << "1. The first parameter indicates the number of cohorts, which should be a positive integer" << endl;
    cout << "2. Any folders in the result save path will not be created, please make sure they exist" << endl;
    cout << "3. Please put the options at the end after result_save_path" << endl << endl;

    cout << "Options:" << endl;
    cout << "-extract: file path of lists of all SNPs included for analysis" << endl;
    cout << "-fixedSNP: file path of lists of all fixed candidate SNPs" << endl;
    cout << "-colinear: colinear_threshold (default: 0.9)" << endl;
    cout << "-R2: R2_incremental_threshold (default: -1, no threshold)" << endl;
    cout << "-R2back: R2_incremental_threshold_backwards (default: -1, no threshold for backward selection)" << endl;
    cout << "-iter_num: total iteration number (default: 10000)" << endl;
    cout << "-window (Mb): SNP position window +/- (default: 10 [which means +/-10Mb], set to -1 for no window" << endl;
    cout << "--no_fast_inv: use normal matrix inverse (default: use fast inverse)" << endl << endl;
    cout << "--no_MDISA: do not run MDISA after MACOJO" << endl;
    
    cout << "Examples:" << endl;
    cout << "./macojo 2 cojoFile1_path PLINK1_path cojoFile2_path PLINK2_path result_save_path -colinear 0.9" << endl;
    cout << "./macojo 1 cojoFile_path PLINK_path result_save_path -fixedSNP user_file -window 10" << endl << endl;
    exit(-1);
}
