cmake_minimum_required(VERSION 3.0)

project(dna-sequencing-alt-negative-errors)

add_executable(heuristic main.cpp ga_solver.cpp individual.cpp xml_reader.cpp)

target_link_libraries(heuristic PUBLIC tinyxml2)

if (CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(heuristic PRIVATE "-Wall" "-Wextra")
endif()
