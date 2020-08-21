
cmake_minimum_required(VERSION 3.10.0)
project(swe CXX)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
Please make a build directory and run CMake from there. Remove CMakeCache.txt.")
endif()

# default to release if build type is empty
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# we need c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenBLAS REQUIRED NO_DEFAULT_PATH)
if (NOT ${OpenBLAS_FOUND})
  message(FATAL_ERROR "OpenBLAS not found")
endif()

find_package(OpenMP REQUIRED)

set(TRILINOS_LIBRARIES kokkosalgorithms
  kokkoscontainers
  kokkoscore
  teuchoskokkoscomm
  teuchoskokkoscompat
  teuchosremainder
  teuchosnumerics
  teuchoscomm
  teuchosparameterlist
  teuchosparser
  teuchoscore
  tpetraext
  tpetrainout
  tpetra
  tpetraclassiclinalg
  tpetraclassicnodeapi
  tpetraclassic
  kokkoskernels
  # repeat to solve issue we have on linux
  kokkosalgorithms
  teuchosparameterlist)

link_directories(${TRILINOS_LIB_DIR})
include_directories(${EIGEN_INCLUDE_DIR} ${TRILINOS_INC_DIR} ${PRESSIO_INCLUDE_DIR})
link_libraries(OpenMP::OpenMP_CXX ${OpenBLAS_LIBRARIES} ${TRILINOS_LIBRARIES})

# fom and rom executables
set(exeName sweFom)
add_executable(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/main_fom.cc)
set(exeName sweRom)
add_executable(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/main_rom.cc)

# python scripts
configure_file(./python_scripts/make_hyper_basis.py
  ${CMAKE_CURRENT_BINARY_DIR}/make_hyper_basis.py COPYONLY)
configure_file(./python_scripts/postProcess.py
  ${CMAKE_CURRENT_BINARY_DIR}/postProcess.py COPYONLY)
