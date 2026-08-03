cmake_minimum_required(VERSION 3.5...3.27)

project( SandiaDecayCascadeCensus )

# Database-wide invariant census. Registered per delivered XML variant, matching
# the xml_file_tests pattern -- one database per invocation keeps peak memory to a
# single rapidxml DOM rather than three.
add_executable( sd_invariant_census sd_invariant_census.cpp )
target_link_libraries( sd_invariant_census SandiaDecay )
target_include_directories( sd_invariant_census PRIVATE ${RAPID_XML_DIR} )
set_property( TARGET sd_invariant_census PROPERTY CXX_STANDARD 11 )

set( _sd_xml_dir "${CMAKE_CURRENT_SOURCE_DIR}/../.." )

add_test( NAME sandia_decay_invariants_full COMMAND sd_invariant_census
  "${_sd_xml_dir}/sandia.decay.xml" --assert --tag full )
add_test( NAME sandia_decay_invariants_min COMMAND sd_invariant_census
  "${_sd_xml_dir}/sandia.decay.min.xml" --assert --tag min )
add_test( NAME sandia_decay_invariants_nocoinc COMMAND sd_invariant_census
  "${_sd_xml_dir}/sandia.decay.nocoinc.min.xml" --assert --tag nocoinc )
