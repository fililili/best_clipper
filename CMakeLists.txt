cmake_minimum_required(VERSION 3.29)
project(best_clipper)

if(NOT DEFINED CPM_SOURCE_CACHE)
  set(CPM_SOURCE_CACHE "$ENV{HOME}/.cache/cpm" CACHE PATH "CPM cache directory")
endif()
include(cmake/CPM.cmake)

CPMAddPackage(
    NAME boost
    GITHUB_REPOSITORY boostorg/boost
    GIT_TAG boost-1.90.0
    OPTIONS
        "BOOST_ENABLE_CMAKE ON"
        "BOOST_INCLUDE_LIBRARIES geometry"
)
get_property(targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "Available targets: ${targets}")

if(TARGET Boost::geometry)
    message(STATUS "Boost::geometry target found")
else()
    message(FATAL_ERROR "Boost::geometry target not found!")
endif()

add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
add_subdirectory(benchmarks)