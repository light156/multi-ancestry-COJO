# Manc-COJO (Multi-ancestry Conditional and Joint Analysis)

**Manc-COJO** is a tool for multi-ancestry conditional and joint analysis of GWAS summary statistics.


## Download

You can either clone the repository:

```bash
git clone https://github.com/light156/multi-ancestry-COJO.git
````

Or simply download the `manc_cojo` binary to your server.

## Usage & Example Commands

```bash
<program_path> <cohort_num> <sumstat_file1> <PLINK_file1> [<sumstat_file2> <PLINK_file2> ...] <output_file> [options]
```
  
For two cohorts, your command should be like
```bash
./manc_cojo 2 GWAS_Cohort1_path LD_reference_Cohort1_path GWAS_Cohort2_path LD_reference_Cohort2_path Output_directory_and_name -colinear 0.9
```

For one cohort, your command should be like
```bash
./manc_cojo 1 GWAS_path LD_reference_path Output_directory_and_name -window 10 --LD -extract Included_SNP_path
```

## Notes

1. The first parameter indicates the number of cohorts, which should be a positive integer.
2. Any folders in the result save path will not be created; please ensure they already exist.
3. Please place all options at the end, after `Output_directory_and_name`.
4. GWAS inputs should be in standard GCTA format: `SNP A1 A2 freq b se N`.
5. LD reference files should be in PLINK format: `.bim`, `.fam`, `.bed` (or `.bim`, `.fam`, `.ld` for using LD inputs, indicated by --LD).

## Options

| Option            | Description                                                                   | Default             |
| ----------------- | ----------------------------------------------------------------------------- | ------------------- |
| `-extract`        | File path for list of SNPs to include in the analysis.                        |                     |
| `-fixedSNP`       | File path for list of fixed candidate SNPs (not removable during iterations). |                     |
| `-colinear`       | Colinearity threshold.                                                        | `0.9`               |
| `-R2`             | R² threshold for forward selection.                                           | `-1` (no threshold) |
| `-R2back`         | R² threshold for backward selection.                                          | `-1` (no threshold) |
| `-window`         | LD window size in Mb. Use `-1` to disable windowing.                          | `10` (±10Mb)        |
| `-freq`           | frequency threshold to exclude rare SNPs.                                     | `0.01`              |
| `--freq_mode_and` | Only keep SNPs that reach frequency threshold in sumstat of all cohorts.      |                     |
|                   | By default, keep SNPs that reaches threshold in at least one cohort.          |                     |
| `--no_MDISA`      | Skip running MDISA after Manc-COJO.                                           |                     |
| `--LD`            | Read PLINK .ld files instead of PLINK .bed files.                             |                     |
| `-iter_num`       | Maximum number of iterations.                                                 | `10000`             |

## A Mock Example: Three Ancestries with HapMap3 SNPs and 20Mb Window; skip running MDISA

```bash
./manc_cojo 3 \
/users/your_name/cohort1_sumstat.linear_gcta_format \
/users/your_name/cohort1_LD_ref \
/users/your_name/cohort2_sumstat.linear_gcta_format \
/users/your_name/cohort2_LD_ref \
/users/your_name/cohort3_sumstat.linear_gcta_format \
/users/your_name/cohort3_LD_ref \
/users/your_name/output_filename \
-colinear 0.9 \
-extract HapMap3.SNPlist \
-window 20 \
--no_MDISA
```

---

Please contact Yong (yong.wang@stats.ox.ac.uk) for software related enquries, or Mark (xiaotong.wang@psych.ox.ac.uk) for algorithm related questions
