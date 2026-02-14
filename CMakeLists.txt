cmake_minimum_required(VERSION 3.29)
project(best_clipper)

if(MSVC)
    add_compile_options("/utf-8")
endif()

set(CPM_SOURCE_CACHE "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH "CPM local packages directory")
set(CPM_DEBUG 1)
include(cmake/CPM.cmake)

CPMAddPackage(
    NAME GraphBLAS
    GITHUB_REPOSITORY DrTimothyAldenDavis/GraphBLAS
    GIT_TAG v10.3.1
    OPTIONS
        "GRAPHBLAS_BUILD_STATIC_LIBS ON"
        "BUILD_SHARED_LIBS OFF"
        "BUILD_STATIC_LIBS ON"
        "GRAPHBLAS_USE_CUDA ON"
)

if(TARGET GraphBLAS_static)
    message(STATUS "GraphBLAS target found")
else()
    message(FATAL_ERROR "GraphBLAS_static target not found!")
endif()

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

add_executable(best_clipper_basic test.cpp)
target_compile_features(best_clipper_basic PRIVATE cxx_std_20)

target_link_libraries(best_clipper_basic PRIVATE GraphBLAS_static Boost::geometry Boost::container)