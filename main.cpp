#include "macojo.h"


int main(int argc, char** argv) 
{
    // double tStart = clock();
    double tStart = omp_get_wtime();

    MACOJO macojo;
    LOGGER << setprecision(12);

    if (argc >= 7 && strcmp(argv[1], "--d") == 0) 
    {   
        LOGGER.open(string(argv[6])+".log");
        
        macojo.initialize_hyperparameters(argc-7, argv+7);
        macojo.read_files_two_cohorts(argv[2], argv[3], argv[4], argv[5]);
        
        LOGGER.i(0, "Loop started");
        macojo.main_loop(argv[6]);
    } 
    else if (argc >= 5 && strcmp(argv[1], "--s") == 0) 
    {   
        LOGGER.open(string(argv[4])+".log");

        macojo.initialize_hyperparameters(argc-5, argv+5);        
        macojo.read_files_one_cohort(argv[2], argv[3]);

        LOGGER.i(0, "Loop started");
        macojo.MDISA(macojo.c1);
        macojo.save_results_DISA(macojo.c1, string(argv[4])+".MDISA.jma.cojo");
    }
    else 
        macojo.show_tips_and_exit();
        
    // LOGGER << "Total running time: " << fixed << setprecision(2) << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "Total running time: " << fixed << setprecision(2) << omp_get_wtime() - tStart << " seconds" << endl;
    LOGGER.close();

    return 0;
}