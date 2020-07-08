cmake_minimum_required(VERSION 3.0)

project(dna-sequencing-alt-negative-errors)

add_subdirectory(tinyxml2-8.0.0)

add_executable(sbh sbh.cpp xml_reader.cpp)

target_link_libraries(sbh PUBLIC tinyxml2)
