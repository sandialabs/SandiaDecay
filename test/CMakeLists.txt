cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project( SandiaDecayTests )

add_executable( xml_file_tests xml_file_tests.cpp )
target_link_libraries( xml_file_tests SandiaDecay )


