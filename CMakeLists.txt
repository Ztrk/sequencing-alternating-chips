# CMAKE_BUILD_TYPE: Release, Debug, Relwithdebinfo, ...

cmake_minimum_required(VERSION 3.0)

project(dna-sequencing-alt-negative-errors)

add_subdirectory(tinyxml2-8.0.0)
add_subdirectory(exact)
add_subdirectory(heuristic)
