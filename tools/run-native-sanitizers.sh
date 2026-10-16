#!/usr/bin/env bash

set -euo pipefail

mode="${1:-}"
case "${mode}" in
  asan|ubsan) ;;
  *)
    echo "usage: tools/run-native-sanitizers.sh asan|ubsan" >&2
    exit 2
    ;;
esac

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/neurocluster-${mode}.XXXXXX")"
lib_dir="${work_dir}/library"
makevars="${work_dir}/Makevars"
base_makevars="${R_MAKEVARS_USER:-${HOME}/.R/Makevars}"
mkdir -p "${lib_dir}"
trap 'rm -rf "${work_dir}"' EXIT

common_flags="-O1 -g -fno-omit-frame-pointer -fno-optimize-sibling-calls"
if [[ "${mode}" == "asan" ]]; then
  sanitizer_flags="-fsanitize=address"
else
  sanitizer_flags="-fsanitize=undefined -fno-sanitize-recover=all"
fi

{
  if [[ -f "${base_makevars}" ]]; then
    printf 'include %s\n' "${base_makevars}"
  fi
  printf 'CXXFLAGS=%s %s\n' "${common_flags}" "${sanitizer_flags}"
  printf 'CXX11FLAGS=%s %s\n' "${common_flags}" "${sanitizer_flags}"
  printf 'CXX14FLAGS=%s %s\n' "${common_flags}" "${sanitizer_flags}"
  printf 'CXX17FLAGS=%s %s\n' "${common_flags}" "${sanitizer_flags}"
  printf 'CXX20FLAGS=%s %s\n' "${common_flags}" "${sanitizer_flags}"
  printf 'LDFLAGS += %s\n' "${sanitizer_flags}"
} > "${makevars}"

export R_MAKEVARS_USER="${makevars}"
export R_LIBS_USER="${lib_dir}${R_LIBS_USER:+:${R_LIBS_USER}}"

if [[ "${mode}" == "asan" ]]; then
  export ASAN_OPTIONS="abort_on_error=1:detect_leaks=0:strict_string_checks=1"
  cxx_config="$(R CMD config CXX17)"
  cxx_bin="${cxx_config%% *}"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    asan_runtime="$("${cxx_bin}" -print-file-name=libclang_rt.asan_osx_dynamic.dylib)"
    if [[ ! -f "${asan_runtime}" ]]; then
      echo "unable to locate the AddressSanitizer runtime via ${cxx_bin}" >&2
      exit 1
    fi
    export DYLD_INSERT_LIBRARIES="${asan_runtime}${DYLD_INSERT_LIBRARIES:+:${DYLD_INSERT_LIBRARIES}}"
  else
    # Do not export LD_PRELOAD here. It would interpose ASan into every child
    # of R CMD INSTALL, including GNU make: with strict_string_checks=1, make's
    # own strndup() reads are reported as heap-buffer-overflows and make aborts
    # (status 134) before any package source is compiled (issue #2). Only the
    # uninstrumented R binary that loads the instrumented package needs the
    # runtime preloaded; the smoke loop below scopes it to that process.
    asan_runtime="$("${cxx_bin}" -print-file-name=libasan.so)"
    if [[ ! -f "${asan_runtime}" ]]; then
      echo "unable to locate the AddressSanitizer runtime via ${cxx_bin}" >&2
      exit 1
    fi
    asan_preload="${asan_runtime}${LD_PRELOAD:+:${LD_PRELOAD}}"
  fi
else
  export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"
fi

# The explicit smoke cases below are the load test. --no-test-load is needed
# on macOS because R CMD INSTALL launches its load probe through a protected
# system shell that strips DYLD_INSERT_LIBRARIES before ASan can interpose.
R CMD INSTALL --preclean --clean --no-multiarch --no-test-load \
  --library="${lib_dir}" "${repo_root}"

r_home="$(R RHOME)"

# Run R's real executable with the ASan runtime preloaded into that process
# only. `R CMD` supplies R's environment (R_HOME, library paths) without the
# preload reaching its own shell helpers.
run_asan_r() {
  R CMD env LD_PRELOAD="${asan_preload}" "${r_home}/bin/exec/R" \
    --vanilla --slave "$@"
}

if [[ "${mode}" == "asan" && "$(uname -s)" != "Darwin" ]]; then
  # Fail loudly rather than pass a smoke matrix that ran without ASan.
  run_asan_r -e 'maps <- readLines("/proc/self/maps"); if (!any(grepl("libasan", maps, fixed = TRUE))) stop("ASan runtime is not loaded in the smoke R process"); cat("asan runtime: loaded\n")'
fi

smoke_cases=(
  slic4d_core
  flash3d_supervoxels_cpp
  fused_assignment
  fused_assignment_parallel_binned
  compute_masked_distances_cpp
  find_1nn_subgraph_cpp
  find_connected_components_cpp
  build_grid_adjacency_cpp
  rena_rnn_coarse_cpp
  ward_on_supervoxels_cpp
  slice_msf_runwise
  slice_fuse_consensus
)

for smoke_case in "${smoke_cases[@]}"; do
  if [[ "${mode}" == "asan" && "$(uname -s)" == "Darwin" ]]; then
    # The Rscript launcher passes through a protected shell on macOS, which
    # strips DYLD_INSERT_LIBRARIES. Invoke R's real executable directly.
    R_HOME="${r_home}" NEUROCLUSTER_USE_INSTALLED=true \
      "${r_home}/bin/exec/R" --vanilla --slave \
        --file="${repo_root}/tools/run-native-subprocess-smoke.R" --args \
        --root="${repo_root}" --case="${smoke_case}"
  elif [[ "${mode}" == "asan" ]]; then
    NEUROCLUSTER_USE_INSTALLED=true \
      run_asan_r --file="${repo_root}/tools/run-native-subprocess-smoke.R" \
        --args --root="${repo_root}" --case="${smoke_case}"
  else
    NEUROCLUSTER_USE_INSTALLED=true \
      Rscript "${repo_root}/tools/run-native-subprocess-smoke.R" \
        --root="${repo_root}" --case="${smoke_case}"
  fi
done

echo "${mode} native smoke: ok (${#smoke_cases[@]} cases)"
