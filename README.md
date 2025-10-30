# Manc-COJO (Multi-ancestry Conditional and Joint Analysis)

**Manc-COJO** is a tool for multi-ancestry conditional and joint analysis (COJO) of GWAS summary statistics.
Note that our program can also perform single-ancestry COJO and reproduce the result of [original GCTA COJO](https://yanglab.westlake.edu.cn/software/gcta/#COJO), but runs much faster.

For example, for HDL trait on ~6,500,000 SNPs and ~76,000 individuals, the running time per chromosome for our program using 1 thread, and GCTA using 5 threads, is as follows. 
![time_comparison_HDL.png](bin_macOS_win/time_comparison_HDL.png)

You can clone this repo or directly download the executable `manc_cojo` for immediate use on 64-bit Intel Linux.  
If it doesn’t run immediately, make sure it has execution permission:

```bash
chmod +x manc_cojo
```

If it still fails to run, your system may not be compatible with the precompiled binary. 
In that case, please follow the build steps (which is quite simple) below or contact us for assistance.

### Verify installation

After building or downloading the executable, you can confirm that it runs correctly by checking its usage information:

```bash
./manc_cojo --help
```

If the program is installed properly, it will print the list of available options and a brief description of each.


> Our program mainly targets Linux servers, but we also include ready-to-use executables for macOS and Windows in the `bin_macOS_win` folder. They were compiled and tested on macOS 15.3.2 and Windows 11. 

> If you run into compatibility issues on any system, feel free to reach out and we are very happy to help.

## Usage & Example Commands

The usage is basically the same as original GCTA COJO, while `--bfile` and `--cojo-file` are extended to multiple cohorts. Nevertheless, there are two major differences:
1. We do not have `--cojo-slct` flag, since the program is designed for performing COJO only.
2. Our model is slightly different from original GCTA version. If you want to **reproduce original GCTA results or use GCTA-COJO criteria on multiple cohorts**, please use flag `--gcta`.

For example, if you want to obtain the **same** result as original GCTA-COJO for a single cohort, your command should be like

```bash
./manc_cojo
--bfile LD_ref_path \
--cojo-file GWAS_sumstat_path \
--out Output_directory_and_name \
--gcta
```

For multiple cohorts with the model in our paper, your command should be like

```bash
./manc_cojo
--bfile LD_ref_path1 LD_ref_path2 ... LD_ref_pathN \
--cojo-file GWAS_sumstat_path1 GWAS_sumstat_path2 ... GWAS_sumstat_pathN \
--out Output_directory_and_name 
```

For calculating joint effects of given SNPs on a single cohort using original GCTA-COJO model

```bash
./manc_cojo \
--bfile myfolder/LD_ref
--cojofile myfolder/sumstat \
--out output_path/output_filename \
--diff-freq 1 \
--cojo-collinear 0.99 \
--extract given_SNP_path \
--cojo-joint --gcta
```

(`--cojo-collinear 0.99` and `--diff-freq 1` for including all given SNPs)


## Notes
1. Please make sure the paths are paired in `--bfile` and `--cojo-file`.
2. For using genotype data, please provide `.bed` `.bim` `.fam` files.
3. For using LD matrix, please provide `.bim` `.ld` files, and add an `--LD` flag in the options. 
   To exclude SNPs with too high MAF difference, please provide `.frq` files as well.
4. The folders in the output file path (which is after `--out`) must exist.
   (You can easily create them by 'mkdir -p <folder_path>')


Please refer to the descriptions for all supported options and flags below.
### Original GCTA Options/Flags 
These options and flags are functionally identical to those in the original GCTA. Please refer to [https://yanglab.westlake.edu.cn/software/gcta/#COJO](https://yanglab.westlake.edu.cn/software/gcta/#COJO) for detailed definitions.

| Option             | Description                                                  | Default          |
| ------------------ | ------------------------------------------------------------ | ---------------- |
| `--bfile`          | PLINK binary file prefix for each cohort                     | `Required`       |
| `--cojo-file`      | GWAS summary statistics file for each cohort                 | `Required`       |
| `--out`            | Output file path prefix                                      | `Required`       |
| `--cojo-wind`      | SNP position window in Kb (`-1` disables windowing)          | `10000` (±1e7 )  |
| `--cojo-p`         | Significance threshold for SNP selection                     | `5e-8`           |
| `--cojo-collinear` | Colinearity threshold for SNP inclusion (`0–0.999`)          | `0.9`            |
| `--maf`            | Minor allele frequency threshold (`1e-5-0.5`)                | `0.01`           |
| `--geno`           | Missingness threshold (`0-1`)                                | `1` (none)       |
| `--diff-freq`      | Frequency diff threshold between sumstat and PLINK (`0-1`)   | `0.2`            |
| `--extract`        | File path for list of SNPs to include in the analysis        |                  |

| Flag               |                                                                 |                  
| ------------------ | --------------------------------------------------------------- |               
| `--cojo-joint`     | Output joint results for given SNPs and exit<br>Only valid when `--extract` is used|

### Multi-ancestry COJO Options/Flags

| Option       | Description                                                       | Default             |
| ------------ | ----------------------------------------------------------------- | ------------------- |
| `--fixed`    | File path for fixed candidate SNPs (non-removable in selection)   |                     |
| `--R2`       | R² threshold for forward selection                                | `-1` (none)         |
| `--R2back`   | R² threshold for backward selection                               | `-1` (none)         |

| Flag               | Description                                                        |
| ------------------ | ------------------------------------------------------------------ | 
| `--freq-mode-and`  | Only keep SNPs that reach MAF threshold in sumstat of all cohorts<br>By default, keep SNPs that reach threshold in at least one cohort |
| `--MDISA`          | Run MDISA after Manc-COJO<br> By default, only run COJO selection on multiple cohorts and exit |

### Main logic Options/Flags

| Option             | Description                                                    | Default       |
| ------------------ | -------------------------------------------------------------- | ------------- |
| `--iter`           | Maximum number of iterations                                   | `10000`       |
| `--thread-num`     | Number of thread to use (One thread is actually fast enough)   | `1`           | 

| Flag               | Description                                                        |
| ------------------ | ------------------------------------------------------------------ | 
| `--gcta`           | Use original GCTA-COJO selection model<br>By default, our model is used|  
| `--LD`             | Read PLINK .ld files instead of PLINK .bed files<br>By default, find .bed files to read|
| `--remove-NA`      | Remove NA in genotype data for correlation calculation, only work for PLINK.bed files<br>By default, use mean imputation on NA genotypes   | 

---
## Third-Party Libraries
This project includes or depends on some third-party open source libraries as follows.

Throughout the program, all matrix computations are based on Eigen 
- [Eigen 3.4.1](https://eigen.tuxfamily.org): external/Eigen

The code for parsing command line options is modified from
- [CLI11](https://github.com/CLIUtils/CLI11): external/CLI11.hpp

The code for logging is modified from
- [GCTA](https://github.com/jianyangqt/gcta/blob/master/include/Logger.h): external/LOGGER.h, external/LOGGER.cpp

---

## Building from Source

If you’d like the fastest performance or plan to modify the code, please clone the repository and build it yourself from the command line. Below are simple example commands for **Linux** and **macOS**.

> **Requirements**
> - A C++11-compatible compiler (`GCC ≥ 4.8.1` or `Clang ≥ 3.3`)
> - Optional: OpenMP for multithreading (default on Linux; not included with Apple Clang)

### Linux

```bash
git clone https://github.com/light156/multi-ancestry-COJO.git
cd multi-ancestry-COJO

g++ -std=c++11 -O3 -march=native -DNDEBUG -fopenmp -pthread \
    -I external/Eigen -I external -I data -I include \
    external/Logger.cpp data/Geno.cpp src/*.cpp \
    -o manc_cojo
```

### macOS
macOS’s default compiler (Apple Clang) does not include OpenMP support. Although you can install LLVM via Homebrew for full OpenMP functionality, using one thread is already fast enough as shown above.
So it is fine to just build a single-threaded version for simplicity:

```bash
git clone https://github.com/light156/multi-ancestry-COJO.git
cd multi-ancestry-COJO

clang++ -std=c++11 -O3 -march=native -DNDEBUG \
    -I external/Eigen -I external -I data -I include \
    external/Logger.cpp data/Geno.cpp src/*.cpp \
    -o manc_cojo
```

### Windows
Since large computing clusters primarily run on Linux, we do not specifically target Windows.
However, if there is a real need to run on Windows, you can install [mingw-w64](https://www.msys2.org/) by following the official guide.
After successful installation, you can compile the source code with GCC using the same commands as on Linux to obtain the executable file.

> Precompiled binaries for Linux, macOS, and Windows are all available in our GitHub repository.

---
Please contact Yong (yong.wang@stats.ox.ac.uk) for software-related enquries and bug reports, or Mark (xiaotong.wang@psych.ox.ac.uk) for algorithm-related questions.
