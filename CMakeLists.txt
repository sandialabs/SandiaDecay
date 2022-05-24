cmake_minimum_required(VERSION 3.1...3.20)

project(SandiaDecay)


set( RAPID_XML_DIR "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty" )
#find_path( RAPID_XML_DIR rapidxml/rapidxml.hpp ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty )

set( sources ${CMAKE_CURRENT_SOURCE_DIR}/SandiaDecay.cpp )
set( headers ${CMAKE_CURRENT_SOURCE_DIR}/SandiaDecay.h )

set( OTHER_SUPPORT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.md )

add_library( SandiaDecay STATIC ${sources} ${headers} ${OTHER_SUPPORT_FILES} )

target_include_directories( SandiaDecay PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                            PRIVATE ${RAPID_XML_DIR} )

#build tools and examples only if building as a standalone project
get_directory_property(hasParent PARENT_DIRECTORY)
if( NOT hasParent  )
  add_subdirectory( tools )
  add_subdirectory( examples )
  add_subdirectory( test )

#  add_definitions( -DALLOW_TRANSITION_MODIFICATION )
#  add_subdirectory( wrappers/c )
endif( NOT hasParent )


