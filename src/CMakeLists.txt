add_library(best_clipper_basic INTERFACE)

target_include_directories(best_clipper_basic INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(best_clipper_basic INTERFACE cxx_std_20)

target_link_libraries(best_clipper_basic INTERFACE GraphBLAS_static Boost::geometry Boost::container)