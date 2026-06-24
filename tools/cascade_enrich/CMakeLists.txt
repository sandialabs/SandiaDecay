cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project( SandiaDecayCascadeEnrich CXX )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

find_path( RAPID_XML_DIR rapidxml/rapidxml.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../../3rdparty
)

# Core G4-parse + angular-correlation + structure-match library (no rapidxml /
# SandiaDecay deps).
add_library( cascade_enrich STATIC
  G4LevelData.cpp
  G4RadDecayData.cpp
  Wigner.cpp
  AngularCorrelation.cpp
  StructureMatch.cpp
)
target_include_directories( cascade_enrich PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

# The XML-injection driver (needs rapidxml).
add_executable( enrich_cascade_data enrich_cascade_data.cpp )
target_link_libraries( enrich_cascade_data cascade_enrich )
target_include_directories( enrich_cascade_data PUBLIC ${RAPID_XML_DIR} )

# Unit tests.
enable_testing()
add_executable( test_g4parsers test_g4parsers.cpp )
target_link_libraries( test_g4parsers cascade_enrich )
add_test( NAME cascade_enrich_parsers COMMAND test_g4parsers )

add_executable( test_angular test_angular.cpp )
target_link_libraries( test_angular cascade_enrich )
add_test( NAME cascade_enrich_angular COMMAND test_angular )

add_executable( test_structure test_structure.cpp )
target_link_libraries( test_structure cascade_enrich )
add_test( NAME cascade_enrich_structure COMMAND test_structure )

# Match-diagnostic tool (audits SD<->G4 gamma matching; needs rapidxml).
add_executable( cascade_match_diag cascade_match_diag.cpp )
target_link_libraries( cascade_match_diag cascade_enrich )
target_include_directories( cascade_match_diag PUBLIC ${RAPID_XML_DIR} )
