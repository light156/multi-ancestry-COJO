#include "macojo.h"


int main(int argc, char** argv) 
{   
    LOGGER << setprecision(12);

    MACOJO macojo;

    // double tStart = clock();
    double tStart = omp_get_wtime();

    if (argc < 2 || atoi(argv[1]) < 1) {
        macojo.show_tips();
        LOGGER.e(0, "Please specify the number of cohorts");
    }
    
    for (int i = 0; i < atoi(argv[1]); i++)
        macojo.cohorts.push_back(Cohort()); 

    int param_num = atoi(argv[1]) * 2 + 3; // cohort num, cojo files, PLINK files, output file 

    if (argc >= param_num) {   
        string savename = argv[param_num-1];
        LOGGER.open(savename+".log");
        
        macojo.read_user_hyperparameters(argc-param_num, argv+param_num);
        macojo.read_cojo_PLINK_files(argv+2, atoi(argv[1]));
        macojo.entry_function(savename);
    } else {
        macojo.show_tips();
        LOGGER.e(0, "Not enough parameters");
    }
        
    // LOGGER << "Total running time: " << fixed << setprecision(2) << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    LOGGER << "Total running time: " << fixed << setprecision(2) << omp_get_wtime() - tStart << " seconds" << endl;
    LOGGER.close();

    return 0;
}