#include "tcojo.h"


void TCOJO::save_results_main_loop(string filepath) 
{   
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    string headers[25] = {"SNP", "A1", "A2", 
        "freq.x", "b.x", "se.x", "p.x", "N.x", "bJ.x", "seJ.x", "zJ.x", "pJ.x", 
        "freq.y", "b.y", "se.y", "p.y", "N.y", "bJ.y", "seJ.y", "zJ.y", "pJ.y", 
        "bJ.ma" , "seJ.ma", "zJ.ma", "pJ.ma"};
    
    for (int i = 0; i < 24; i++) 
        jmaCOJO << headers[i] << "\t";
        
    jmaCOJO << headers[24] << "\n";

    int index = 0;
    map<string, int> SNP_order_pair;
 
    // Inserting element in pair vector to keep track of previous indexes
    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, index++) 
        SNP_order_pair.insert(make_pair(final_commonSNP[*iter], index));
    
    ArrayXd temp_row;
    ArrayXd bJx = c1.output_b, bJy = c2.output_b;
    ArrayXd seJx = sqrt(c1.output_se2), seJy = sqrt(c2.output_se2);
    ArrayXd zJx = bJx / seJx, zJy = bJy / seJy;
    ArrayXd pJx = erfc(abs(zJx)/sqrt(2)), pJy = erfc(abs(zJy)/sqrt(2));
    
    ArrayXd bJma = sumstat_merge.col(0);
    ArrayXd seJma = sqrt(sumstat_merge.col(1));
    ArrayXd zJma = bJma / seJma;
    ArrayXd pJma = erfc(abs(zJma)/sqrt(2));
    
    jmaCOJO.precision(12);

    for (auto iter = SNP_order_pair.begin(); iter != SNP_order_pair.end(); iter++) {
        index = iter->second;
        jmaCOJO << iter->first << "\t" << A1_ref[final_commonSNP_index[candidate_SNP[index]]] << "\t" \
            << A2_ref[final_commonSNP_index[candidate_SNP[index]]] << "\t";

        // cols 0:b, 1:se2, 2:p, 3:freq, 4:N
        temp_row = c1.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t";
        jmaCOJO << bJx(index) << "\t" << seJx(index) << "\t" << zJx(index) << "\t" << pJx(index) << "\t";
            
        temp_row = c2.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t";
        jmaCOJO << bJy(index) << "\t" << seJy(index) << "\t" << zJy(index) << "\t" << pJy(index) << "\t";
        
        jmaCOJO << bJma(index) << "\t" << seJma(index) << "\t" << zJma(index) << "\t" << pJma(index) << "\n";
    }
    jmaCOJO.clear();
    jmaCOJO.close();

    map<string, int>().swap(SNP_order_pair);
    LOGGER.i(0, "Results saved into [" + filepath + "]\n", "Finished");
}


void TCOJO::save_results_DISA(Cohort &c, string filepath) 
{   
    ofstream jmaCOJO(filepath.c_str());
    if (!jmaCOJO) LOGGER.e(0, "cannot open the file [" + filepath + "] to write.");
    
    string headers[12] = {"SNP", "A1", "A2", "freq", "b", "se", "p", "N", "bJ", "seJ", "zJ", "pJ"};
        
    for (int i = 0; i < 11; i++) 
        jmaCOJO << headers[i] << "\t";
        
    jmaCOJO << headers[11] << "\n";

    int index = 0;
    map<string, int> SNP_order_pair;
 
    // Inserting element in pair vector to keep track of previous indexes
    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, index++) 
        SNP_order_pair.insert(make_pair(final_commonSNP[*iter], index));
    
    ArrayXd temp_row;
    ArrayXd bJ = c.output_b;
    ArrayXd seJ = sqrt(c.output_se2);
    ArrayXd zJ = bJ / seJ;
    ArrayXd pJ = erfc(abs(zJ)/sqrt(2));

    jmaCOJO.precision(12);

    for (auto iter = SNP_order_pair.begin(); iter != SNP_order_pair.end(); iter++) {
        index = iter->second;
        jmaCOJO << iter->first << "\t" << A1_ref[final_commonSNP_index[candidate_SNP[index]]] << "\t" \
            << A2_ref[final_commonSNP_index[candidate_SNP[index]]] << "\t";

        // cols 0:b, 1:se2, 2:p, 3:freq, 4:N
        temp_row = c.sumstat_candidate.row(index);
        jmaCOJO << temp_row(3) << "\t" << temp_row(0) << "\t" << sqrt(temp_row(1)) << "\t" 
            << temp_row(2) << "\t" << temp_row(4) << "\t";

        jmaCOJO << bJ(index) << "\t" << seJ(index) << "\t" << zJ(index) << "\t" << pJ(index) << "\n";
    }
    jmaCOJO.clear();
    jmaCOJO.close();

    map<string, int>().swap(SNP_order_pair);
    LOGGER.i(0, "Results saved into [" + filepath + "]\n", "Finished");
}


void TCOJO::initialize_hyperparameters(int argc, char** argv) 
{
    int temp_num = 0;

    while (temp_num < argc) {
        if (strcmp(argv[temp_num], "-colinear") == 0 && temp_num+1 < argc) {
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
            window_size = LONG_MAX;
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "-window") == 0 && temp_num+1 < argc && atoi(argv[temp_num+1]) > 0) {
            window_size = atoi(argv[temp_num+1]);
            temp_num += 2;
        } else if (strcmp(argv[temp_num], "--no_fast_inv") == 0) {
            if_fast_inv = false;
            temp_num += 1;
        } else 
            show_tips_and_exit();
    }

    colinear_threshold_sqrt = sqrt(colinear_threshold);
    iter_colinear_threshold = 1 / (1 - colinear_threshold);

    LOGGER << "Threshold: 5e-8" << endl;
    LOGGER << "Colinearity threshold: " << colinear_threshold << endl;
    LOGGER << "R2 incremental threshold: " << R2_incremental_threshold << endl;
    LOGGER << "R2 incremental threshold backwards: " << R2_incremental_threshold_backwards << endl;
    LOGGER << "Maximal iteration number: " << max_iter_num << endl;
    LOGGER << "SNP position window (+/-): " << window_size << endl;
    if (if_fast_inv)
        LOGGER << "Use fast matrix inversion" << endl;
    else
        LOGGER << "Use normal matrix inversion" << endl;

    LOGGER << "--------------------------------" << endl;
}


void TCOJO::show_tips_and_exit() 
{
    cout << endl << "Program usage: Joint meta analysis for two cohorts (--d) or a single cohort (--s)" << endl;
    cout << "Usage 1: program_path --d cojoFile1_path PLINK1_path cojoFile2_path PLINK2_path result_save_path (+other options)" << endl;
    cout << "Usage 2: program_path --s cojoFile_path PLINK_path result_save_path (+other options)" << endl << endl;

    cout << "TIPS:" << endl;
    cout << "1. The first parameter must be --d or --s indicating the number of cohorts" << endl;
    cout << "2. Any folders in the result save path will not be created, please make sure they exist" << endl;
    cout << "3. Please put the options at the end. Possible options include (default value in brackets):" << endl;
    cout << "-colinear: colinear_threshold (0.9)" << endl;
    cout << "-R2: R2_incremental_threshold (0)" << endl;
    cout << "-R2back: R2_incremental_threshold_backwards (-0.5)" << endl;
    cout << "-iter_num: total iteration number (10000)" << endl;
    cout << "-window: SNP position window +/- (1500000), set to -1 for no window" << endl;
    cout << "--no_fast_inv: use normal matrix inverse" << endl << endl;
    
    cout << "For example, both commands below will give the same result, as they use default option values:" << endl;
    cout << "./tcojo --s cojoFile_path PLINK_path result_save_path -colinear 0.9 -R2back -0.5 -iter_num 10000" << endl;
    cout << "./tcojo --s cojoFile_path PLINK_path result_save_path -R2 0 " << endl << endl;
    exit(-1);
}
