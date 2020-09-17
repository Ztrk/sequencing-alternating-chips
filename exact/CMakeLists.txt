cmake_minimum_required(VERSION 3.0)

project(dna-sequencing-alt-negative-errors)

add_subdirectory(tinyxml2-8.0.0)

add_executable(main main.cpp ga_solver.cpp xml_reader.cpp DNAPath.cpp overlapGraph.cpp)

target_link_libraries(main PUBLIC tinyxml2)
