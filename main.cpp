#include "trans_ancestry_cojo.h"


void show_tips_and_exit() 
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


int main(int argc, char** argv) {
    
    clock_t tStart = clock();

    TransAncestryCOJO tacojo;

    LOGGER << setprecision(12);

    if (argc >= 7 && strcmp(argv[1], "--d") == 0) 
    {   
        LOGGER.open(string(argv[6])+".log");
        
        bool flag = tacojo.initialize_hyperparameters(argc-7, argv+7);
        if (!flag) show_tips_and_exit();
            
        tacojo.read_files_two_cohorts(argv[2], argv[3], argv[4], argv[5]);

        LOGGER.i(0, "Loop started");
        tacojo.main_loop(argv[6]);
    } 
    else if (argc >= 5 && strcmp(argv[1], "--s") == 0) 
    {   
        LOGGER.open(string(argv[4])+".log");

        bool flag = tacojo.initialize_hyperparameters(argc-5, argv+5);
        if (!flag) show_tips_and_exit();

        tacojo.read_files_one_cohort(argv[2], argv[3]);

        LOGGER.i(0, "Loop started");
        tacojo.MDISA(tacojo.c1);
        tacojo.save_results_DISA(tacojo.c1, string(argv[4])+".MDISA.jma.cojo");
    }
    else 
        show_tips_and_exit();
        
    LOGGER << "Total running time: " << fixed << setprecision(2) << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER.close();

    return 0;
}