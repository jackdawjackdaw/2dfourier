cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
project (2dfourier)

find_package(GSL REQUIRED)
message(STATUS "GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_INSTALL_PATH: ${CMAKE_INSTALL_PREFIX}")

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
MESSAGE(STATUS "include dirs: ${CMAKE_include_directories}")
SET(libs ${libs} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
MESSAGE(STATUS "libs: ${libs}")
add_subdirectory(src)
