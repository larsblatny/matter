
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)


if (USE_VDB)

    add_executable(tests
      tests.cpp
      tests_vdb.cpp
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
      ../src/simulation/nonlocal_projection.cpp
			../src/simulation/nonlocal_p2g.cpp
			../src/simulation/nonlocal_g2p.cpp
    )

	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

	find_package(OpenVDB 3.1 QUIET)
	include_directories(${OPENVDB_INCLUDE_DIRS} SYSTEM)

	find_package(TBB REQUIRED)
	include_directories(${TBB_INCLUDE_DIRS} SYSTEM)

	find_package(Boost REQUIRED)
	include_directories(${Boost_INCLUDE_DIRS} SYSTEM)

  target_link_libraries (tests GTest::gtest_main
                      OpenMP::OpenMP_CXX
                      Eigen3::Eigen
                      ${OPENVDB_LIBRARIES}
                      ${OpenEXR_LIBRARIES}
                      ${TBB_LIBRARIES}
                      ${Boost_LIBRARIES}
                  )

else (USE_VDB)

    add_executable(tests
      tests.cpp
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

	target_link_libraries (tests GTest::gtest_main
    OpenMP::OpenMP_CXX
    Eigen3::Eigen)

endif (USE_VDB)

target_compile_definitions(tests PRIVATE INCLUDE_DIR="${CMAKE_SOURCE_DIR}/levelsets")
target_compile_definitions(tests PRIVATE SRC_DIR="${CMAKE_SOURCE_DIR}/src")

get_target_property(DEFS mpm COMPILE_DEFINITIONS)
message(STATUS "Compile Definitions: ${DEFS}")

include(GoogleTest)
gtest_discover_tests(tests)
