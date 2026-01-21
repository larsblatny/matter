find_package (Eigen3 3.3 REQUIRED NO_MODULE)
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
			   simulation/nonlocal_projection.cpp
			   simulation/nonlocal_p2g.cpp
			   simulation/nonlocal_g2p.cpp
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