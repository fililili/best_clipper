add_library(best_clipper_core
    best_clipper.cpp
    chain_builder.cpp
    half_chain_graph.cpp
    connect_half_chains.cpp
    output_builder.cpp
    graph_helper.cpp
    filter_half_chains.cpp
)

target_include_directories(best_clipper_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(best_clipper_core PUBLIC cxx_std_20)

target_link_libraries(best_clipper_core PUBLIC Boost::geometry Boost::container Boost::multiprecision)
