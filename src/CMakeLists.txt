# Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

# cmake_policy(SET CMP0167 OLD)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_executable(mpm
			   mpm.cpp
			   simulation/simulation.cpp
			   simulation/save_data.cpp
			   simulation/remesh.cpp
			   simulation/deformation_update.cpp
			   simulation/plasticity.cpp
			   simulation/p2g.cpp
			   simulation/g2p.cpp
			   simulation/explicit_euler_update.cpp
			   simulation/boundary_collision.cpp
			   simulation/update_dt.cpp
			   simulation/position_update.cpp
			   simulation/musl.cpp
			   simulation/pbc.cpp
			   )

# Enable LTO for mpm only (not GoogleTest) to avoid CMP0069 warnings from third-party targets
if(ipo_supported)
    set_property(TARGET mpm PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    message(STATUS "IPO/LTO enabled for mpm")
endif()

# SPECIFIC RELEASE OPTIMIZATIONS (safer fast-math subset — preserves NaN/Inf detection)
target_compile_options(mpm PRIVATE
    $<$<CONFIG:Release>:-fno-math-errno>
    $<$<CONFIG:Release>:-fno-trapping-math>
    $<$<CONFIG:Release>:-fno-signed-zeros>
    $<$<CONFIG:Release>:-fassociative-math>
    $<$<CONFIG:Release>:-freciprocal-math>
    $<$<CONFIG:Release>:-funroll-loops>
    $<$<CONFIG:Release>:-fomit-frame-pointer>
)

# EIGEN-SPECIFIC OPTIMIZATIONS
target_compile_definitions(mpm PRIVATE
    EIGEN_NO_DEBUG                    # Disable debug checks in Eigen for speed
    EIGEN_UNROLLING_LIMIT=256         # Unroll small fixed-size matrix loops (raised from 128)
    EIGEN_FAST_MATH=1                 # Use fast math approximations where possible
    EIGEN_VECTORIZE_FMA               # Enable fused multiply-add instructions
)

if (USE_VDB)

	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

	find_package(OpenVDB 3.1 QUIET)
	include_directories(${OPENVDB_INCLUDE_DIRS} SYSTEM)

	find_package(TBB REQUIRED)
	include_directories(${TBB_INCLUDE_DIRS} SYSTEM)

	find_package(Boost REQUIRED)
	include_directories(${Boost_INCLUDE_DIRS} SYSTEM)

	target_link_libraries (mpm
						OpenMP::OpenMP_CXX
						Eigen3::Eigen
						${OPENVDB_LIBRARIES}
						${OpenEXR_LIBRARIES}
						${TBB_LIBRARIES}
					${Boost_LIBRARIES}
					)

else (USE_VDB)

	target_link_libraries (mpm OpenMP::OpenMP_CXX Eigen3::Eigen)

endif (USE_VDB)

target_compile_definitions(mpm PRIVATE INCLUDE_DIR="${CMAKE_SOURCE_DIR}/levelsets")
target_compile_definitions(mpm PRIVATE SRC_DIR="${CMAKE_SOURCE_DIR}/src")

get_target_property(DEFS mpm COMPILE_DEFINITIONS)
message(STATUS "Compile Definitions: ${DEFS}")