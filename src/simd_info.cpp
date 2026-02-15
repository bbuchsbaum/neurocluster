// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>

#include <cstdlib>
#include <string>

using namespace Rcpp;

namespace {

static std::string compiler_id() {
#if defined(__clang__)
  return std::string("clang ") + std::to_string(__clang_major__) + "." +
         std::to_string(__clang_minor__) + "." + std::to_string(__clang_patchlevel__);
#elif defined(__GNUC__)
  return std::string("gcc ") + std::to_string(__GNUC__) + "." +
         std::to_string(__GNUC_MINOR__) + "." + std::to_string(__GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER)
  return std::string("msvc ") + std::to_string(_MSC_VER);
#else
  return "unknown";
#endif
}

static std::string arch_id() {
#if defined(__aarch64__) || defined(_M_ARM64)
  return "aarch64";
#elif defined(__arm__) || defined(_M_ARM)
  return "arm";
#elif defined(__x86_64__) || defined(_M_X64)
  return "x86_64";
#elif defined(__i386__) || defined(_M_IX86)
  return "x86";
#else
  return "unknown";
#endif
}

static const char *dot_i8_kernel_name() {
  // This is a stable API surface: once we add true ISA-specific kernels, we
  // will switch this string based on build flags / dispatch.
#if defined(__aarch64__) && defined(__ARM_FEATURE_DOTPROD)
  const char *val = std::getenv("NEUROCLUSTER_FORCE_PORTABLE_DOT_I8");
  const bool forced = (val && val[0] != '\0' && val[0] != '0');
  return forced ? "portable_forced" : "neon_dotprod";
#endif
#if defined(NEUROCLUSTER_DOT_I8_AVX512VNNI)
  return "avx512vnni";
#elif defined(NEUROCLUSTER_DOT_I8_AVX2)
  return "avx2";
#elif defined(NEUROCLUSTER_DOT_I8_NEON_DOTPROD)
  return "neon_dotprod";
#else
#ifdef _OPENMP
  return "portable_omp_simd";
#else
  return "portable";
#endif
#endif
}

static List compile_features() {
  return List::create(
    _["openmp"] =
#ifdef _OPENMP
      true,
#else
      false,
#endif
    _["avx2"] =
#if defined(__AVX2__)
      true,
#else
      false,
#endif
    _["avx512f"] =
#if defined(__AVX512F__)
      true,
#else
      false,
#endif
    _["avx512bw"] =
#if defined(__AVX512BW__)
      true,
#else
      false,
#endif
    _["avx512vnni"] =
#if defined(__AVX512VNNI__)
      true,
#else
      false,
#endif
    _["sse4_1"] =
#if defined(__SSE4_1__)
      true,
#else
      false,
#endif
    _["neon"] =
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
      true,
#else
      false,
#endif
    _["arm_dotprod"] =
#if defined(__ARM_FEATURE_DOTPROD)
      true
#else
      false
#endif
  );
}

static List runtime_features() {
  // Keep runtime detection conservative to ensure compilation everywhere.
  // "Real SIMD kernels" ultimately means ISA-specific instructions are present
  // in the binary; use the disassembly checker script for that.
  List out = List::create(
    _["available"] = false
  );

#if (defined(__clang__) || defined(__GNUC__)) && (defined(__x86_64__) || defined(__i386__))
  // __builtin_cpu_supports requires string literals; keep to widely-supported ones.
  out["available"] = true;
  out["sse4_1"] = static_cast<bool>(__builtin_cpu_supports("sse4.1"));
  out["avx2"] = static_cast<bool>(__builtin_cpu_supports("avx2"));
#endif

  return out;
}

} // namespace

// [[Rcpp::export]]
List simd_info_cpp() {
  return List::create(
    _["compiler"] = compiler_id(),
    _["arch"] = arch_id(),
    _["compile"] = compile_features(),
    _["runtime"] = runtime_features(),
    _["dot_i8_kernel"] = dot_i8_kernel_name(),
    _["notes"] =
      "To verify ISA-specific int8 SIMD dot kernels, disassemble the shared "
      "library and search for dot-product instructions. See "
      "inst/benchmarks/check_simd_int8.R."
  );
}
