cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project( SandiaDecayTools )

find_path( RAPID_XML_DIR rapidxml/rapidxml.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty
  ${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty
)

add_executable( minimize_xml_file minimize_xml_file.cpp )
target_link_libraries( minimize_xml_file SandiaDecay )
target_include_directories( minimize_xml_file PUBLIC ${RAPID_XML_DIR} )

add_executable( decay_database_benchmark decay_database_benchmark.cpp )
target_link_libraries( decay_database_benchmark SandiaDecay )
target_include_directories( decay_database_benchmark PUBLIC ${RAPID_XML_DIR} )

