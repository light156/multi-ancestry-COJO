#include "doctest.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "../include/utils.hpp"
#include "../include/macojo.h"

static std::string make_temp_path(const std::string& name) {
    // TMPDIR is set by test/run_tests.sh; fallback to current dir if not set.
    const char* tmpdir = std::getenv("TMPDIR");
    std::string base = (tmpdir && *tmpdir) ? tmpdir : ".";
    return base + "/" + name;
}

TEST_CASE("skim_fam reads FID/IID pairs") {
    const std::string path = make_temp_path("fam_test.fam");
    {
        std::ofstream out(path.c_str());
        out << "F1 I1 0 0 1 -9\n";
        out << "F2 I2 0 0 2 -9\n";
    }

    std::vector<std::string> ids;
    skim_fam(path, ids);
    CHECK(ids.size() == 2);
    CHECK(ids[0] == "F1:I1");
    CHECK(ids[1] == "F2:I2");

    std::remove(path.c_str());
}

TEST_CASE("skim_bim filters by chromosome and extracts SNPs") {
    const std::string path = make_temp_path("bim_test.bim");
    {
        std::ofstream out(path.c_str());
        out << "1 rs1 0 100 A G\n";
        out << "2 rs2 0 200 C T\n";
        out << "1 rs3 0 300 G A\n";
    }

    std::vector<std::string> snps_chr1;
    skim_bim(path, 1, snps_chr1);
    CHECK(snps_chr1.size() == 2);
    CHECK(snps_chr1[0] == "rs1");
    CHECK(snps_chr1[1] == "rs3");

    std::vector<std::string> snps_chr2;
    skim_bim(path, 2, snps_chr2);
    CHECK(snps_chr2.size() == 1);
    CHECK(snps_chr2[0] == "rs2");

    std::remove(path.c_str());
}

TEST_CASE("skim_SNP reads SNP list with header and column index") {
    const std::string path = make_temp_path("snp_test.txt");
    {
        std::ofstream out(path.c_str());
        out << "id snp other\n";
        out << "1 rs10 x\n";
        out << "2 rs20 y\n";
    }

    std::vector<std::string> snps;
    skim_SNP(std::vector<std::string>{path, "header", "2"}, snps);
    CHECK(snps.size() == 2);
    CHECK(snps[0] == "rs10");
    CHECK(snps[1] == "rs20");

    std::remove(path.c_str());
}

// Stub to satisfy linker when read_files.cpp references output_bad_SNP.
void MACOJO::output_bad_SNP(string) {}
