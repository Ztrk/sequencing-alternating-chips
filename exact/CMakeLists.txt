# CMAKE_BUILD_TYPE: Release, Debug, ...

cmake_minimum_required(VERSION 3.0)

project(dna-sequencing-alt-negative-errors)

add_executable(exact main.cpp exact_solver.cpp xml_reader.cpp DNAPath.cpp overlapGraph.cpp)

target_link_libraries(exact PUBLIC tinyxml2)
 
if (CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(exact PRIVATE "-Wall" "-Wextra")
endif()
