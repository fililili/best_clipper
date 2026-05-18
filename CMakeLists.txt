cmake_minimum_required(VERSION 3.29)
project(best_clipper)

if(MSVC)
    add_compile_options("/utf-8")
endif()

set(CPM_SOURCE_CACHE "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH "CPM local packages directory")
set(CPM_DEBUG 1)
include(cmake/CPM.cmake)

CPMAddPackage(
    NAME boost
    GITHUB_REPOSITORY boostorg/boost
    GIT_TAG boost-1.90.0
)
get_property(targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "Available targets: ${targets}")

if(TARGET Boost::geometry)
    message(STATUS "Boost::geometry target found")
else()
    message(FATAL_ERROR "Boost::geometry target not found!")
endif()

if(TARGET Boost::container)
    message(STATUS "Boost::container target found")
else()
    message(FATAL_ERROR "Boost::container target not found!")
endif()

# oneTBB (required by Embree)
CPMAddPackage(
    NAME oneTBB
    GITHUB_REPOSITORY oneapi-src/oneTBB
    GIT_TAG v2021.13.0
    OPTIONS "TBB_TEST OFF" "TBB_EXAMPLES OFF" "TBB_STRICT OFF" "TBBMALLOC_BUILD OFF"
)
if(TARGET tbb AND NOT TARGET TBB::tbb)
    add_library(TBB::tbb ALIAS tbb)
endif()
if(TARGET tbbmalloc AND NOT TARGET TBB::tbbmalloc)
    add_library(TBB::tbbmalloc ALIAS tbbmalloc)
endif()
set(TBB_FOUND TRUE CACHE BOOL "TBB found via CPM" FORCE)

# Intel Embree - high-performance ray tracing
CPMAddPackage(
    NAME embree
    GITHUB_REPOSITORY embree/embree
    GIT_TAG v4.3.3
    OPTIONS
    "EMBREE_ISPC_SUPPORT OFF"
    "EMBREE_TUTORIALS OFF"
    "EMBREE_TESTING OFF"
    "EMBREE_STATIC_LIB ON"
)

add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
add_subdirectory(benchmarks)
add_subdirectory(debug)