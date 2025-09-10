# Manc-COJO (Multi-ancestry Conditional and Joint Analysis)

**Manc-COJO** is a tool for multi-ancestry conditional and joint analysis of GWAS summary statistics.


## Download

You can either clone the repository:

```bash
git clone https://github.com/light156/multi-ancestry-COJO.git
````

Or simply download the `manc-cojo` binary to your server.

## Example Command

```bash
./macojo 2 GWAS_Cohort1_path LD_reference_Cohort1_path GWAS_Cohort2_path LD_reference_Cohort2_path Output_directory_and_name -colinear 0.9
```

## Notes

1. The first parameter indicates the number of cohorts, which should be a positive integer.
2. Any folders in the result save path will not be created; please ensure they already exist.
3. Please place the options at the end, after `Output_directory_and_name`.
4. GWAS inputs should be in standard GCTA format: `SNP A1 A2 freq b se N`.
5. LD reference files should be in PLINK format: `.bim`, `.fam`, `.bed`.

## Options

| Option          | Description                                                                   | Default             |
| --------------- | ----------------------------------------------------------------------------- | ------------------- |
| `-extract`      | File path for list of SNPs to include in the analysis.                        | —                   |
| `-fixedSNP`     | File path for list of fixed candidate SNPs (not removable during iterations). | —                   |
| `-colinear`     | Colinearity threshold.                                                        | `0.9`               |
| `-R2`           | R² threshold for forward selection.                                           | `-1` (no threshold) |
| `-R2back`       | R² threshold for backward selection.                                          | `-1` (no threshold) |
| `-iter_num`     | Maximum number of iterations.                                                 | `10000`             |
| `-window`       | LD window size in Mb. Use `-1` to disable windowing.                          | `10`                |
| `--no_fast_inv` | Disable the fast matrix inversion algorithm.                                  | —                   |
| `--no_MDISA`    | Skip running MDISA after Manc-COJO.                                           | —                   |

## Example: Three Ancestries with HapMap3 SNPs and 20Mb Window; turn off fast matrix inversion algorithm; skip running MDISA

```bash
./macojo 3 \
GWAS_Cohort1_path \
LD_reference_Cohort1_path \
GWAS_Cohort2_path \
LD_reference_Cohort2_path \
GWAS_Cohort3_path \
LD_reference_Cohort3_path \
Output_directory_and_name \
-colinear 0.9 \
-extract HapMap3.SNPlist \
-window 20 \
--no_fast_inv \
--no_MDISA
```

---

Please contact Yong (yong.wang@stats.ox.ac.uk) for software related enquries, or Mark (xiaotong.wang@psych.ox.ac.uk) for algorithm related questions
