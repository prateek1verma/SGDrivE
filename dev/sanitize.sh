#!/usr/bin/env bash
# Local AddressSanitizer + UndefinedBehaviorSanitizer run of the C++ tests on
# macOS (Apple clang), a quick alternative to the manual `sanitizers` GitHub
# workflow. Builds the package with sanitizers into a temporary library and
# runs the C++-related test files with the ASan runtime preloaded. The R
# binary is started directly because macOS drops DYLD_* variables when the
# `R` shell script is used.
#
# Run from the package root:  bash dev/sanitize.sh
set -euo pipefail
[ -f DESCRIPTION ] || { echo "run from the package root"; exit 1; }

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir "$tmp/lib"
cat > "$tmp/Makevars" <<'MK'
CXX17FLAGS = -g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=undefined
LDFLAGS = -fsanitize=address,undefined
MK
# build from a tarball so that no sanitizer objects are left in src/
R CMD build --no-build-vignettes --no-manual . >"$tmp/build.log" 2>&1 || { cat "$tmp/build.log"; exit 1; }
mv ScrewDrivE_*.tar.gz "$tmp/"
R_MAKEVARS_USER="$tmp/Makevars" R CMD INSTALL --no-test-load -l "$tmp/lib" "$tmp"/ScrewDrivE_*.tar.gz \
  >"$tmp/install.log" 2>&1 || { cat "$tmp/install.log"; exit 1; }

asan="$(clang++ -print-file-name=libclang_rt.asan_osx_dynamic.dylib)"
r_home="$(R RHOME)"
cat > "$tmp/run.R" <<RS
.libPaths(c("$tmp/lib", .libPaths()))
res <- testthat::test_dir("tests/testthat", package = "ScrewDrivE", load_package = "installed",
  filter = "rates-cpp|dopri5|model-view|simulate|backend-equivalence|regression|stress|cb-",
  reporter = "summary", stop_on_failure = TRUE)
RS
NOT_CRAN=true ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=print_stacktrace=1 R_HOME="$r_home" \
  DYLD_INSERT_LIBRARIES="$asan" "$r_home/bin/exec/R" --vanilla --quiet -f "$tmp/run.R"
echo "No sanitizer errors."
