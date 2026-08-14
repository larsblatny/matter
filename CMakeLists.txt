# Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

cmake_minimum_required(VERSION 3.16)
SET (CMAKE_CXX_STANDARD 17)
project(matter)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ===== OPTIMIZATION FLAGS FOR RELEASE BUILD =====
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Add aggressive optimization flags
    if(MSVC)
        # MSVC compiler flags
        add_compile_options(/O2 /GL /arch:AVX2)
        message(STATUS "Compile options: /O2 /GL /arch:AVX2")
        string(APPEND CMAKE_EXE_LINKER_FLAGS " /LTCG")
    else()
        # GCC/Clang compiler flags
        add_compile_options(-O3 -march=native -mtune=native -ftree-vectorize -fvariable-expansion-in-unroller)
        message(STATUS "Compile options: -O3 -march=native -mtune=native -ftree-vectorize -fvariable-expansion-in-unroller")
        # Suppress false positive allocation warnings from OpenVDB LTO analysis
        string(APPEND CMAKE_CXX_FLAGS " -Wno-alloc-size-larger-than")
    endif()
endif()

# Check LTO/IPO support — result used below to enable it on both the mpm executable
# and the matter Python module
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output)
if(NOT ipo_supported)
    message(STATUS "IPO/LTO not supported: ${ipo_output}")
endif()
# =================================================

# Compile settings shared between the mpm executable (src/) and the matter Python
# module (python/), so simulation code is compiled identically either way — single
# source of truth, referenced by both src/CMakeLists.txt and python/CMakeLists.txt.
set(MATTER_EIGEN_DEFS
    EIGEN_NO_DEBUG                    # Disable debug checks in Eigen for speed
    EIGEN_UNROLLING_LIMIT=256         # Unroll small fixed-size matrix loops (raised from 128)
    EIGEN_FAST_MATH=1                 # Use fast math approximations where possible
    EIGEN_VECTORIZE_FMA               # Enable fused multiply-add instructions
)
set(MATTER_RELEASE_FLAGS
    # SPECIFIC RELEASE OPTIMIZATIONS (safer fast-math subset — preserves NaN/Inf detection)
    $<$<CONFIG:Release>:-fno-math-errno>
    $<$<CONFIG:Release>:-fno-trapping-math>
    $<$<CONFIG:Release>:-fno-signed-zeros>
    $<$<CONFIG:Release>:-fassociative-math>
    $<$<CONFIG:Release>:-freciprocal-math>
    $<$<CONFIG:Release>:-funroll-loops>
    $<$<CONFIG:Release>:-fomit-frame-pointer>
)

option (USE_VDB "Use OpenVDB part of the code" ON)
message(STATUS "Use of OpenVDB: ${USE_VDB}")

option (BUILD_PYTHON_BINDINGS "Build pybind11 Python bindings for Simulation" ON)
message(STATUS "Build Python bindings: ${BUILD_PYTHON_BINDINGS}")

enable_testing()

# Custom target to clean CMake cache (like make clean for cache)
add_custom_target(clean-cache
    COMMAND ${CMAKE_COMMAND} -E remove -f CMakeCache.txt CMakeFiles/cmake.check_cache
    COMMENT "Removing CMake cache files"
)

add_subdirectory(src)
add_subdirectory(tests)

if (BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

#### Download header files from external Github repos ####

# Create deps folder
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/deps/)

# Download tinyply header file
file(DOWNLOAD
    https://raw.githubusercontent.com/ddiakopoulos/tinyply/refs/heads/master/source/tinyply.h
    ${CMAKE_SOURCE_DIR}/deps/tinyply.h)

# Download tph_possion v0.3
file(DOWNLOAD
    https://github.com/thinks/tph_poisson/archive/refs/tags/v0.3.tar.gz
    ${CMAKE_SOURCE_DIR}/deps/tph_poisson.tar.gz)
# Extract zipped folder
file(ARCHIVE_EXTRACT INPUT 
    ${CMAKE_SOURCE_DIR}/deps/tph_poisson.tar.gz
    DESTINATION ${CMAKE_SOURCE_DIR}/deps/)
# Delete the zipped folder
file(REMOVE ${CMAKE_SOURCE_DIR}/deps/tph_poisson.tar.gz)