cmake_minimum_required(VERSION 3.13)
project(ULB)

set(CMAKE_CXX_STANDARD 14)

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
            "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

option(USE_OMP "use OpenMP" ON)
option(USE_MPI "use MPI" ON)
option(BUILD_PYINT "build python interface" ON)
option(BUILD_TEST "build tests" ON)
option(BUILD_TETGEN "build tetgen library" ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(cpp)
if(BUILD_PYINT)
    add_subdirectory(swig)
endif()
if(BUILD_TEST)
    add_subdirectory(test)
endif()
if(BUILD_TETGEN)
    add_subdirectory(lib/tetgen)
    add_subdirectory(lib/triangle)
    include_directories(lib)
endif()