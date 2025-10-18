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
<program_path> <cohort_num> <sumstat_file1> <PLINK_file1> [<sumstat_file2> <PLINK_file2> ...] --out <output_file> [options]
```
  
For two cohorts, your command should be like
```bash
./manc_cojo 2 \
GWAS_Cohort1_path LD_ref_Cohort1_path \
GWAS_Cohort2_path LD_ref_Cohort2_path \
--out Output_directory_and_name
```

For one cohort, your command should be like
```bash
./manc_cojo 1 GWAS_path LD_ref_path \
--out Output_directory_and_name
```

## Notes

1. The first parameter indicates the number of cohorts, which should be a positive integer.
2. The next 2*N parameters are file paths of GWAS summary statistics and PLINK files for N cohorts, respectively.
3. All options must be placed after the above parameters.
4. The directory specified in the output path (after --out) must already exist.
   (You can create it in advance, for example with 'mkdir -p <folder_path>')
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
|                   | By default, keep SNPs that reach threshold in at least one cohort.            |                     |
| `--no_MDISA`      | Skip running MDISA after Manc-COJO.                                           |                     |
| `--LD`            | Read PLINK .ld files instead of PLINK .bed files.                             |                     |
| `--cojo-joint`    | Only output for provided fixed candidate SNPs and exit.                       |                     |
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
--out /users/your_name/output_filename \
--cojo-collinear 0.99 \
--extract HapMap3.SNPlist \
--window 20
```

---

Please contact Yong (yong.wang@stats.ox.ac.uk) for software related enquries, or Mark (xiaotong.wang@psych.ox.ac.uk) for algorithm related questions
