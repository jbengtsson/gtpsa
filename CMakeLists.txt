cmake_minimum_required(VERSION 3.5...4.0)

project("gtpsa" LANGUAGES C CXX)
set(GTPSA_VERSION 0.0.1)

# -----------------------------------------------------------------------------
# CRITICAL FIX:
# Force use of in-tree gtpsa headers before any system-installed gtpsa
# -----------------------------------------------------------------------------
include_directories(BEFORE
  ${CMAKE_CURRENT_SOURCE_DIR}/c++
)

OPTION (USE_OpenMP "Use OpenMP" ON)
OPTION (USE_pybind11 "Use cmake to build the pybind11 extension" ON)
option(PYTHON_MODULE_DEVELOPMENT_INSTALL "Install python binary module to source directory" ON)

find_package(BLAS)
find_package(LAPACK)
find_package(Armadillo REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
  message(STATUS "Default to Release build")
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(gtpsa_DEBUG_BUILD ON)
else()
  set(gtpsa_DEBUG_BUILD OFF)
endif()

IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    message(STATUS "Open MP enabled using flags CFLAGS=${OpenMP_C_FLAGS} CXXFLAGS=${OpenMP_CXX_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ")
  ENDIF()
ENDIF()

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

enable_testing()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

add_subdirectory(c)
add_subdirectory(c++)

IF(USE_pybind11)
  add_subdirectory(python)
endif()

add_subdirectory(examples)
add_subdirectory(tests)
