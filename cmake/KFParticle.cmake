set(KFParticle_INSTALL_DIR ${EXTERNAL_INSTALL_DIR})
set(KFParticle_INCLUDE_DIR ${KFParticle_INSTALL_DIR}/include)
set(KFParticle_LIBRARY_DIR ${KFParticle_INSTALL_DIR}/lib)

# KFParticle's own CMakeLists.txt hardcodes `-ffast-math` into its `Release` flags
# (`-O3 -ftree-vectorize -ffast-math -DNODEBUG`). `-ffast-math` lets the compiler
# reassociate/vectorize floating-point math, which changes KFParticle's transport /
# chi2 numerical results. That is enough to shift ~15% of event-mixing background
# candidates across the cut thresholds versus the reference (AFKF) build
# (see test/cluster_xcheck/FINDINGS.md).
#
# To keep PFSimple numerically consistent with the reference build, we compile the
# bundled KFParticle with a build type that does NOT add `-ffast-math`. Opt back in
# explicitly with `-DPFSimple_KFParticle_FAST_MATH=ON` (not recommended).
set(PFSimple_KFParticle_FAST_MATH OFF CACHE BOOL
    "Build the bundled KFParticle with -ffast-math (changes FP results; not recommended)")

# Pick the build type to forward to the KFParticle sub-build. `Release` (and an
# empty type, which the top-level defaults to Release) is where KFParticle injects
# `-ffast-math`, so use `RelWithDebInfo` (i.e. `-O2 -g -DNDEBUG`, no `-ffast-math`).
set(KFParticle_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
if(NOT PFSimple_KFParticle_FAST_MATH)
    if(KFParticle_BUILD_TYPE STREQUAL "" OR KFParticle_BUILD_TYPE STREQUAL "Release")
        set(KFParticle_BUILD_TYPE "RelWithDebInfo")
    endif()
endif()

ExternalProject_Add(KFParticle_Ext
        GIT_REPOSITORY  "https://git.cbm.gsi.de/pwg-c2f/analysis/KFParticle.git"
        GIT_TAG         "cmake"
        UPDATE_DISCONNECTED ${UPDATE_DISCONNECTED}
        SOURCE_DIR      "${EXTERNAL_DIR}/KFParticle_src"
        BINARY_DIR      "${EXTERNAL_DIR}/KFParticle_build"
        INSTALL_DIR     "${KFParticle_INSTALL_DIR}"
        CMAKE_ARGS
        "-DEXTERNAL_INSTALL_DIR=${KFParticle_INSTALL_DIR}"
        "-DCMAKE_INSTALL_PREFIX=${KFParticle_INSTALL_DIR}"
        "-DCMAKE_BUILD_TYPE=${KFParticle_BUILD_TYPE}"
        "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
        "-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}"
        "-DFIXTARGET=TRUE"
        )
list(APPEND PROJECT_DEPENDENCIES KFParticle_Ext)
list(APPEND PROJECT_LINK_DIRECTORIES ${KFParticle_LIBRARY_DIR})
list(APPEND PROJECT_INCLUDE_DIRECTORIES ${KFParticle_INCLUDE_DIR})
