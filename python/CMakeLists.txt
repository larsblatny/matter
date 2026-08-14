# Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

include(FetchContent)
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.13.6.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(pybind11)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

pybind11_add_module(matter
               bindings.cpp
               ../src/simulation/simulation.cpp
               ../src/simulation/save_data.cpp
               ../src/simulation/remesh.cpp
               ../src/simulation/deformation_update.cpp
               ../src/simulation/plasticity.cpp
               ../src/simulation/p2g.cpp
               ../src/simulation/g2p.cpp
               ../src/simulation/explicit_euler_update.cpp
               ../src/simulation/boundary_collision.cpp
               ../src/simulation/update_dt.cpp
               ../src/simulation/position_update.cpp
               ../src/simulation/musl.cpp
               ../src/simulation/pbc.cpp
               )

target_link_libraries(matter PRIVATE OpenMP::OpenMP_CXX Eigen3::Eigen)

target_compile_definitions(matter PRIVATE
    EIGEN_NO_DEBUG
    EIGEN_UNROLLING_LIMIT=256
    EIGEN_FAST_MATH=1
    EIGEN_VECTORIZE_FMA
)

if (USE_VDB)

    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

    find_package(OpenVDB 3.1 QUIET)
    include_directories(${OPENVDB_INCLUDE_DIRS} SYSTEM)

    find_package(TBB REQUIRED)
    include_directories(${TBB_INCLUDE_DIRS} SYSTEM)

    find_package(Boost REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS} SYSTEM)

    target_compile_definitions(matter PRIVATE USE_VDB)

    target_link_libraries(matter PRIVATE
                        ${OPENVDB_LIBRARIES}
                        ${OpenEXR_LIBRARIES}
                        ${TBB_LIBRARIES}
                        ${Boost_LIBRARIES}
                    )

endif (USE_VDB)

target_compile_definitions(matter PRIVATE INCLUDE_DIR="${CMAKE_SOURCE_DIR}/levelsets")
target_compile_definitions(matter PRIVATE SRC_DIR="${CMAKE_SOURCE_DIR}/src")
