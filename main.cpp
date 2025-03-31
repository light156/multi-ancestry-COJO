#include "tcojo.h"


int main(int argc, char** argv) 
{
    /*
    #pragma omp parallel for
    for (int i = 0; i < 100; i++) 
    {
        cout << "Thread " << omp_get_thread_num() << endl;
    }
    exit(0);
    */

    double tStart = clock();

    TCOJO tcojo;
    LOGGER << setprecision(12);

    if (argc >= 7 && strcmp(argv[1], "--d") == 0) 
    {   
        LOGGER.open(string(argv[6])+".log");
        
        tcojo.initialize_hyperparameters(argc-7, argv+7);
        tcojo.read_files_two_cohorts(argv[2], argv[3], argv[4], argv[5]);
        
        LOGGER.i(0, "Loop started");
        tcojo.main_loop(argv[6]);
    } 
    else if (argc >= 5 && strcmp(argv[1], "--s") == 0) 
    {   
        LOGGER.open(string(argv[4])+".log");

        tcojo.initialize_hyperparameters(argc-5, argv+5);        
        tcojo.read_files_one_cohort(argv[2], argv[3]);

        LOGGER.i(0, "Loop started");
        tcojo.MDISA(tcojo.c1);
        tcojo.save_results_DISA(tcojo.c1, string(argv[4])+".MDISA.jma.cojo");
    }
    else 
        tcojo.show_tips_and_exit();
        
    LOGGER << "Total running time: " << fixed << setprecision(2) << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER.close();

    return 0;
}