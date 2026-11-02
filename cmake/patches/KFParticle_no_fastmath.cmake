# KFParticle_no_fastmath.cmake
#
# Removes `-ffast-math` from the bundled KFParticle's CMakeLists.txt so the
# KFParticle sub-build does not use it, regardless of the chosen build type.
#
# This is invoked by the KFParticle ExternalProject PATCH_COMMAND when
# `-DPFSimple_KFParticle_NO_FAST_MATH=ON` is set. It edits the freshly-cloned
# KFParticle source before its configure step, so the flag is stripped from
# KFParticle's own flags (which it hardcodes into `CMAKE_CXX_FLAGS_RELEASE`).
#
# Required variable:
#   _KFP_SRC_DIR  - path to the KFParticle source directory.

if(NOT DEFINED _KFP_SRC_DIR OR _KFP_SRC_DIR STREQUAL "")
    message(FATAL_ERROR "KFParticle_no_fastmath: _KFP_SRC_DIR is not set.")
endif()

set(_kfp_cmake "${_KFP_SRC_DIR}/CMakeLists.txt")
if(NOT EXISTS "${_kfp_cmake}")
    message(FATAL_ERROR "KFParticle_no_fastmath: ${_kfp_cmake} not found.")
endif()

file(READ "${_kfp_cmake}" _kfp_content)
string(REGEX REPLACE "-ffast-math[ \t]*" "" _kfp_content "${_kfp_content}")
file(WRITE "${_kfp_cmake}" "${_kfp_content}")

message(STATUS "KFParticle_no_fastmath: removed -ffast-math from ${_kfp_cmake}")
