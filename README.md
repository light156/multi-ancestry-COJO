Welcome to use Manc-COJO (Multi-ancestry Conditional and Joint analysis)

Download the software:
You can either use git clone, or simply download the "manc-cojo" file to your server.

Example Command Line:
./macojo 2 GWAS_Cohort1_path LD_reference_Cohort1_path GWAS_Cohort2_path LD_reference_Cohort2_path output_directory_and_name -colinear 0.9

NOTES:
1. The first parameter indicates the number of cohorts, which should be a positive integer
2. Any folders in the result save path will not be created; please make sure they exist
3. Please put the options at the end, after output_directory_and_name
4. GWAS inputs should be in standard GCTA format (SNP A1 A2 freq b se N)
5. LD reference should be in PLINK format (.bim .fam .bed)

Options:
-extract: file path for lists of all SNPs included for analysis
-fixedSNP: file path for lists of all fixed candidate SNPs that can not be removed during the following iterations
-colinear: colinear_threshold (default: 0.9)
-R2: R2 threshold for forward selection (default: -1, no threshold)
-R2back: R2 threshold for backward selection (default: -1, no threshold for backward selection)
-iter_num: maximum iteration number (default: 10000)
-window: LD window size in Mb unit (default: 10 [which means +/-10Mb], set to -1 if you do not want to include a window)
--no_fast_inv: When inverting the correlation matrix, we use a simplified algorithm, which accelerates the matrix inversion with minimal impact on accuracy. You can turn off this algorithm. (Default: use fast inverse)
--no_MDISA: do not run MDISA after Manc-COJO

Examples command line to run Manc-COJO on three ancestries with only Hapmap3 SNPs, using a 20Mb window:
./macojo 2 \
GWAS_Cohort1_path \
LD_reference_Cohort1_path \
GWAS_Cohort2_path \
LD_reference_Cohort2_path \
GWAS_Cohort3_path \
LD_reference_Cohort3_path \
output_directory_and_name \
-colinear 0.9
-extract HapMap3.SNPlist
-window 20
