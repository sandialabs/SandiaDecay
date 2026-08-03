cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project( SandiaDecayTests )

enable_testing()

add_executable( xml_file_tests xml_file_tests.cpp )
target_link_libraries( xml_file_tests SandiaDecay )
add_test( NAME sandia_decay_xml_full COMMAND xml_file_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.xml" )
add_test( NAME sandia_decay_xml_min COMMAND xml_file_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.min.xml" )
add_test( NAME sandia_decay_xml_nocoinc COMMAND xml_file_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.nocoinc.min.xml" )

add_executable( cascade_detail_tests cascade_detail_tests.cpp )
target_link_libraries( cascade_detail_tests SandiaDecay )
target_include_directories( cascade_detail_tests PRIVATE ${RAPID_XML_DIR} )
target_compile_definitions( cascade_detail_tests PRIVATE
  SANDIA_DECAY_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.." )
add_test( NAME sandia_decay_cascade_detail COMMAND cascade_detail_tests )
