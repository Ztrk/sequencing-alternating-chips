cmake_minimum_required(VERSION 3.0)

project(dna-sequencing-alt-negative-errors)

add_subdirectory(tinyxml2-8.0.0)

add_executable(main main.cpp ga_solver.cpp individual.cpp xml_reader.cpp)

target_link_libraries(main PUBLIC tinyxml2)

if (CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(main PRIVATE "-Wall" "-Wextra")
endif()
