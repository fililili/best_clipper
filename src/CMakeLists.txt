add_library(best_clipper_core INTERFACE)

target_include_directories(best_clipper_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(best_clipper_core INTERFACE cxx_std_20)

target_link_libraries(best_clipper_core INTERFACE Boost::geometry Boost::container Boost::multiprecision)