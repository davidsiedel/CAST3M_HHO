cmake_policy(SET CMP0048 NEW)
include(cmake/modules/mechhano.cmake)
project(mechhcano
        VERSION 0.0.1
        DESCRIPTION "C++ implementation of the Hybrid High Order method"
        HOMEPAGE_URL "https://github.com/davidsiedel/mechhano-0"
        LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)

enable_language (Fortran)


#compiler options
include(cmake/modules/compiler.cmake)
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  include(cmake/modules/CodeCoverage.cmake)
  SETUP_TARGET_FOR_COVERAGE(coverage ctest coverage)
endif(CMAKE_BUILD_TYPE STREQUAL "Coverage")

# use eigen
find_package (Eigen3 REQUIRED NO_MODULE)

# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()

add_subdirectory(src)
add_subdirectory(bindings)

#if (TEST_MODE STREQUAL "test")
#  add_subdirectory(tests)
#endif()
