cmake_minimum_required(VERSION 3.22)

project(AFAC_solver LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()
enable_testing()

include(FetchContent)
FetchContent_Declare(
  test_drive
  GIT_REPOSITORY https://github.com/fortran-lang/test-drive.git
  GIT_TAG v0.6.1
)
FetchContent_MakeAvailable(test_drive)

FetchContent_Declare(
  forbenchmark
  GIT_REPOSITORY https://github.com/gha3mi/forbenchmark.git
  GIT_TAG v0.2.0
)
FetchContent_MakeAvailable(forbenchmark)
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/tests)
set(MAIN_DIR ${CMAKE_SOURCE_DIR}/main)
set(BENCHMARK_DIR ${CMAKE_SOURCE_DIR}/benchmarking)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -cpp -O3 -traceback  -warn all -xHost -qopenmp -coarray=distributed -mcmodel=medium -shared-intel")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fpp -cpp -traceback -warn all -xHost -qopenmp -coarray=distributed -mcmodel=medium -shared-intel")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fpp -cpp -traceback -warn all -xHost -qopenmp -coarray=distributed -mcmodel=medium -shared-intel")

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(main)
add_subdirectory(benchmarking)
