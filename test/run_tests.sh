#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${repo_root}/test"
tmp_dir="${build_dir}/tmp"
test_bin="${build_dir}/run_tests"

mkdir -p "${tmp_dir}"
export TMPDIR="${tmp_dir}"

# Build tests in debug-friendly mode for faster compilation.
g++ -std=c++11 -O0 -g -I "${repo_root}" -I"${repo_root}/test" -I"${repo_root}/include" -I"${repo_root}/external" \
  "${repo_root}/test/test_utils.cpp" \
  "${repo_root}/test/test_read_files.cpp" \
  "${repo_root}/test/test_calc_inv.cpp" \
  "${repo_root}/test/test_geno.cpp" \
  "${repo_root}/src/read_files.cpp" \
  "${repo_root}/src/geno.cpp" -o "${test_bin}"

echo "Built ${test_bin}"
echo "Running tests..."
"${test_bin}"

rm -rf "${tmp_dir}"
rm -f "${test_bin}"
