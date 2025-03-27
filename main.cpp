#include "trans_ancestry_cojo.h"
#include <time.h>


int main(int argc, char** argv) {

    LOGGER << setprecision(12);
    
    clock_t tStart = clock();
    
    /*
    if (argc != 6) {
        cout << "Program needs 5 arguments to run" << endl;
        cout << "Example: ./main cojoFile1_path cojoFile2_path PLINK1_path PLINK2_path result_save_path" << endl;
        return 0;
    }

    string cojoFile1 = argv[1], cojoFile2 = argv[2], PLINK1 = argv[3], PLINK2 = argv[4], save_path = argv[5];
    */

    string cojoFile1 = "Height_AFR_GIANT_QCed_withoutUKB_Chr22.sumstat";
    string cojoFile2 = "Height_EUR_GIANT_QCed_withoutUKB_Chr22.sumstat";
    string PLINK1 = "1000G_Chr22_AFR_QCed";
    string PLINK2 = "1000G_Chr22_EUR_QCed";
    
    // string num = argv[1];
    // cojoFile1 = "Simstat/Chr22_LD_Block_1_Casual_SNPs_20_Scenario_"+num+"_AFR.PHENO1.glm.linear_gcta_format";
    // cojoFile2 = "Simstat/Chr22_LD_Block_1_Casual_SNPs_20_Scenario_"+num+"_EUR1.PHENO1.glm.linear_gcta_format";
    
    TransAncestryCOJO tacojo;
    tacojo.read_files(cojoFile1, cojoFile2, PLINK1, PLINK2);
    tacojo.main_loop("results.test");

    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
}