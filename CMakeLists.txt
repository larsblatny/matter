cmake_minimum_required(VERSION 3.16)
SET (CMAKE_CXX_STANDARD 17)
project(matter)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

option (USE_VDB "Use OpenVDB part of the code" OFF)

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