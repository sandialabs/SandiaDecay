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

# Reviewed A=234 decay data of the U-238 chain (intensities, line ownership, x-rays, coincidences).
add_executable( u238_chain_tests u238_chain_tests.cpp )
target_link_libraries( u238_chain_tests SandiaDecay )
add_test( NAME sandia_decay_u238_chain_full COMMAND u238_chain_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.xml" )
add_test( NAME sandia_decay_u238_chain_min COMMAND u238_chain_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.min.xml" )
add_test( NAME sandia_decay_u238_chain_nocoinc COMMAND u238_chain_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.nocoinc.min.xml" )

# Reviewed low-energy U-235 / U-238 data (gamma energies and intensities, Th K x-ray group totals).
add_executable( u235_lowE_tests u235_lowE_tests.cpp )
target_link_libraries( u235_lowE_tests SandiaDecay )
add_test( NAME sandia_decay_u235_lowE_full COMMAND u235_lowE_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.xml" )
add_test( NAME sandia_decay_u235_lowE_min COMMAND u235_lowE_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.min.xml" )
add_test( NAME sandia_decay_u235_lowE_nocoinc COMMAND u235_lowE_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.nocoinc.min.xml" )

# Th, Pa and U K x-ray structure (decay sets of Th/Pa/U daughters and the Pb/Th/Pa/U fluorescence blocks).
add_executable( xray_k_structure_tests xray_k_structure_tests.cpp )
target_link_libraries( xray_k_structure_tests SandiaDecay )
add_test( NAME sandia_decay_xray_k_structure_full COMMAND xray_k_structure_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.xml" )
add_test( NAME sandia_decay_xray_k_structure_min COMMAND xray_k_structure_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.min.xml" )
add_test( NAME sandia_decay_xray_k_structure_nocoinc COMMAND xray_k_structure_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.nocoinc.min.xml" )

# Individually reviewed decay lines outside the U-235 / U-238 reviews (Sr-89 / Y-89m 909 keV, ...).
add_executable( reviewed_lines_tests reviewed_lines_tests.cpp )
target_link_libraries( reviewed_lines_tests SandiaDecay )
add_test( NAME sandia_decay_reviewed_lines_full COMMAND reviewed_lines_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.xml" )
add_test( NAME sandia_decay_reviewed_lines_min COMMAND reviewed_lines_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.min.xml" )
add_test( NAME sandia_decay_reviewed_lines_nocoinc COMMAND reviewed_lines_tests
  "${CMAKE_CURRENT_SOURCE_DIR}/../sandia.decay.nocoinc.min.xml" )
