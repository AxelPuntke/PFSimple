set(KFParticle_INSTALL_DIR ${EXTERNAL_INSTALL_DIR})
set(KFParticle_INCLUDE_DIR ${KFParticle_INSTALL_DIR}/include)
set(KFParticle_LIBRARY_DIR ${KFParticle_INSTALL_DIR}/lib)

# By default the bundled KFParticle is built exactly as upstream, i.e. with
# whatever flags its own CMakeLists.txt picks for the chosen build type
# (including `-ffast-math` for a `Release` build). Passing
# `-DPFSimple_KFParticle_NO_FAST_MATH=ON` applies a patch that strips
# `-ffast-math` from KFParticle's flags, independent of the build type. This is
# only needed when PFSimple must reproduce the numerical results of a reference
# build that was compiled without `-ffast-math` (e.g. comparing against `afkf`).
# See test/cluster_xcheck/FINDINGS.md.
option(PFSimple_KFParticle_NO_FAST_MATH
    "Patch the bundled KFParticle to remove -ffast-math (default OFF; changes KFParticle numerics)"
    OFF)

set(_pfsimple_kfp_patch "")
if(PFSimple_KFParticle_NO_FAST_MATH)
    set(_pfsimple_kfp_patch
        PATCH_COMMAND ${CMAKE_COMMAND}
            "-D_KFP_SRC_DIR=${EXTERNAL_DIR}/KFParticle_src"
            "-P" "${CMAKE_CURRENT_LIST_DIR}/patches/KFParticle_no_fastmath.cmake")
endif()

ExternalProject_Add(KFParticle_Ext
        GIT_REPOSITORY  "https://git.cbm.gsi.de/pwg-c2f/analysis/KFParticle.git"
        GIT_TAG         "cmake"
        UPDATE_DISCONNECTED ${UPDATE_DISCONNECTED}
        SOURCE_DIR      "${EXTERNAL_DIR}/KFParticle_src"
        BINARY_DIR      "${EXTERNAL_DIR}/KFParticle_build"
        INSTALL_DIR     "${KFParticle_INSTALL_DIR}"
        ${_pfsimple_kfp_patch}
        CMAKE_ARGS
        "-DEXTERNAL_INSTALL_DIR=${KFParticle_INSTALL_DIR}"
        "-DCMAKE_INSTALL_PREFIX=${KFParticle_INSTALL_DIR}"
        "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
        "-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}"
        "-DFIXTARGET=TRUE"
        )
list(APPEND PROJECT_DEPENDENCIES KFParticle_Ext)
list(APPEND PROJECT_LINK_DIRECTORIES ${KFParticle_LIBRARY_DIR})
list(APPEND PROJECT_INCLUDE_DIRECTORIES ${KFParticle_INCLUDE_DIR})
