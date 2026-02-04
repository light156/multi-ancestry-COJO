# Manc-COJO <a href="https://github.com/SMARTbiomed/software-review/issues/3"><img align="right" alt="SMB peer reviewed" src="https://smartbiomed.github.io/software-review/data/smartbiomed-software.svg" /></a>

**Manc-COJO** is a C++ software tool for multi-ancestry conditional and joint analysis (COJO) of GWAS summary statistics. 

GWAS tests for association between a trait and SNPs one at a time, giving marginal SNP effect estimates. 
However, associations detected in GWAS are often not independent because of linkage disequilibrium (LD) between SNPs. 
To address this challenge, [COJO](https://www.nature.com/articles/ng.2213) has been proposed and widely used for single-ancestry analyses, where it identifies independent association signals through iterative conditioning on significant SNPs while jointly modelling their effects to account for LD. 
Building upon COJO, our multi-ancestry extension exploits population-specific LD differences to improve the detection of independent association signals and reduce false positives compared to single-ancestry COJO (and _ad hoc_ adaptations for multi-ancestry use).

<br>

Note that **Manc-COJO** can also perform single-ancestry COJO and reproduce the results of [GCTA-COJO](https://yanglab.westlake.edu.cn/software/gcta/#COJO), while running **substantially faster**.

For example, for the analysis of high-density lipoprotein (HDL) cholesterol in our study, involving approximately 6.5 million SNPs and 76,000 individuals, we compared the per-chromosome runtime of Manc-COJO (using 1 thread) with GCTA-COJO (using 5 threads), as shown below. 

![time_HDL.png](https://light156.github.io/multi-ancestry-COJO-docs/time_HDL.png)

## Installation

You can download our software directly from the links below:

### [⬇️ Download link (Linux)](https://github.com/light156/multi-ancestry-COJO/releases/download/v1.0.0/manc_cojo)

### [⬇️ Download link (Windows)](https://github.com/light156/multi-ancestry-COJO/releases/download/v1.0.0/manc_cojo_win.exe)

### [⬇️ Download link (macOS)](https://github.com/light156/multi-ancestry-COJO/releases/download/v1.0.0/manc_cojo_macOS)

<br>

On **Linux** and **macOS**, you may need to give the file permission to run:

```bash
chmod +x manc_cojo
```

To check that the software is working, run:

```bash
./manc_cojo --help
```

This will print the available options and confirm that it runs correctly.

If you have any problems running the software on your system, please feel free to contact us and we are happy to help.

## Usage

If you are not familiar with COJO, we strongly recommend that you follow [**our tutorial**](https://light156.github.io/multi-ancestry-COJO-docs/tutorial). It includes example data and step-by-step commands.

### Single cohort

The usage is largely consistent with GCTA-COJO. 
In general, you can replace the path to the GCTA executable with `manc_cojo` while keeping all remaining arguments unchanged:

```bash
./manc_cojo \
--bfile path \
--cojo-file GWAS_sumstat_path \
--out Output_path_name \
--cojo-slct
```

### Multiple cohorts

Append multiple file paths to `--bfile` and ``--cojo-file``. Please make sure that the files are provided in the same order for both options.

```bash
./manc_cojo \
--bfile path1 path2 ... pathN \
--cojo-file GWAS_sumstat_path1 GWAS_sumstat_path2 ... GWAS_sumstat_pathN \
--out Output_path_name \
--cojo-slct 
```

<br>

Despite being largely similar, the following behaviours intentionally differ from GCTA-COJO:

1. Our software excludes SNPs whose genotypes are identical across all individuals.
2. In output files, both **A1** and **A2** are reported for each SNP. **A1** corresponds to **refA** in GCTA outputs.  
3. By default, our software does **not** generate `.cma.cojo` and `.ldr.cojo` files, as they can be very large and are not required for most use cases. Use `--output-all` to enable all output files, which will also record unqualified SNPs in the corresponding `.badsnps` files.
4. When collinearity issues arise among user-provided SNPs during conditional analysis (`--cojo-cond`) or joint analysis (`--cojo-joint`), GCTA terminates without output. In contrast, our software iteratively removes problematic SNPs until the issue is resolved. Removed SNPs are recorded in the `.log` file.

<br>

For advanced usage, a complete list of command-line options, and instructions for running on UKB-RAP, please refer to our documentation website: https://light156.github.io/multi-ancestry-COJO-docs

## Citation

If you find our paper or software useful for your research, please consider citing our paper: 

> Multi-ancestry conditional and joint analysis (Manc-COJO) applied to GWAS summary statistics. Xiaotong Wang, Yong Wang, Peter M Visscher, Naomi R Wray, Loic Yengo. bioRxiv 2026.01.30.702783; doi: [https://doi.org/10.64898/2026.01.30.702783](https://doi.org/10.64898/2026.01.30.702783) 

## Support

Please contact Yong (yong.wang@psych.ox.ac.uk) for software-related enquries and bug reports, or Mark (xiaotong.wang@psych.ox.ac.uk) for algorithm-related questions. 
We welcome GitHub issues, including usage feedback and new feature requests, so that discussions are visible to all users.

---

## License and Acknowledgments

This project is released under the **MIT License** (see the `LICENSE` file for details).

This project includes or depends on the following third-party open-source libraries:

- **[Eigen 3.4.1](https://eigen.tuxfamily.org)** – Used for all matrix computations (`external/Eigen`)
- **[CLI11](https://github.com/CLIUtils/CLI11)** – Modified for parsing command-line options (`external/CLI11.hpp`)
- **[fastfloat](https://github.com/fastfloat/fast_float)** - Used for number parsing (`external/fast_float.h`)
- **[doctest](https://github.com/doctest/doctest)** - **Only** for unit tests, **not** linked or used by the main program (`test/doctest.h`)
