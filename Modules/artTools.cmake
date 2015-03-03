# - Developer/User tools for building art-style products
#-----------------------------------------------------------------------
# Copyright 2014 Ben Morgan <Ben.Morgan@warwick.ac.uk>
# Copyright 2014 University of Warwick

# - Include zis only once
if(__ART_TOOLS_INCLUDED)
  return()
endif()
set(__ART_TOOLS_INCLUDED 1)


# - Interfaces
include(CMakeParseArguments)

#-----------------------------------------------------------------------
# Functions to replace "simple_plugin"
# Follow cmake add_<target> type interface, possibly allow for SOURCES
# and LINKS_LIBRARIES style arguments
#
# Art understands three main types of plugin
# "service", "module" and "source".
# Because the search for plugins is based on the library name being in
# the form
#  libAA_bb_Cc_<type>.<ext>
# you have to know what <type> is
# In addition, each specific type of module requires linking to certain
# art libraries.

#-----------------------------------------------------------------------
# function art_add_service(<name> source1 source2 ... sourceN)
#
function(art_add_service _name)
  if(NOT (_name MATCHES ".*_service$"))
    message(FATAL_ERROR "art_add_service: target name must end in '_service'")
  endif()

  add_library(${_name} SHARED ${ARGN})
  target_link_libraries(${_name}
    ${art_IMPORT_NAMESPACE}art_Framework_Services_Registry
    FNALCore::FNALCore
    )
endfunction()


#-----------------------------------------------------------------------
# function art_add_module(<name> source1 source2 ... sourceN)
#
function(art_add_module _name)
  if(NOT (_name MATCHES ".*_module$"))
    message(FATAL_ERROR "art_add_module: target name must end in '_module'")
  endif()

  add_library(${_name} SHARED ${ARGN})
  target_link_libraries(${_name}
    ${art_IMPORT_NAMESPACE}art_Framework_Core
    ${art_IMPORT_NAMESPACE}art_Framework_Principal
    ${art_IMPORT_NAMESPACE}art_Persistency_Common
    ${art_IMPORT_NAMESPACE}art_Persistency_Provenance
    ${art_IMPORT_NAMESPACE}art_Utilities
    FNALCore::FNALCore
    ${ROOT_Core_Library}
    )
endfunction()

#-----------------------------------------------------------------------
# function art_add_source(<name> source1 source2 ... sourceN)
#
function(art_add_source _name)
  if(NOT (_name MATCHES ".*_source$"))
    message(FATAL_ERROR "art_add_source: target name must end in '_source'")
  endif()

  add_library(${_name} SHARED ${ARGN})
  target_link_libraries(${_name}
    ${art_IMPORT_NAMESPACE}art_Framework_Core
    ${art_IMPORT_NAMESPACE}art_Framework_Principal
    ${art_IMPORT_NAMESPACE}art_Framework_IO_Sources
    ${art_IMPORT_NAMESPACE}art_Persistency_Common
    ${art_IMPORT_NAMESPACE}art_Persistency_Provenance
    ${art_IMPORT_NAMESPACE}art_Utilities
    FNALCore::FNALCore
    ${ROOT_Core_Library}
    )
endfunction()


#-----------------------------------------------------------------------
# function art_get_genreflex_flags(<var>)
#          Set value of <var> to minimal genreflex flags for generating
#          dictionaries from Art sources
#function(art_get_genreflex_flags __outputvar)
#  set(${__outputvar}
#    --deep
#    --iocomments
#    --fail_on_warnings
#    --gccxmlopt=--debug
#    --gccxmlopt=--gccxml-compiler
#    --gccxmlopt=${CMAKE_CXX_COMPILER}
#    -D_REENTRANT
#    -DGNU_SOURCE
#    -DGNU_GCC
#    -DPROJECT_NAME="${PROJECT_NAME}"
#    -DPROJECT_VERSION="${art_VERSION}"
#    -D__STRICT_ANSI__
#
#    PARENT_SCOPE
#    )
#endfunction()


#-----------------------------------------------------------------------
# function art_add_dictionary(<name> <selection> <sources>)
#          Build Reflex based dictionary libraries from input headers
#          and selection file
#          A "dictionary" is actually two shared libs <name>_dict and
#          <name>_map
#
# - Inputs: Name, headers, selection file
#function(art_add_dictionary __name __selection)

  # - Selection file may be relative to current dir or absolute
  #  if(IS_ABSOLUTE "${__selection}")
#    set(selection_file "${__selection}")
#  else()
#    set(selection_file "${CMAKE_CURRENT_SOURCE_DIR}/${__selection}")
#  endif()
#
#  if(NOT EXISTS "${selection_file}")
#    message(FATAL_ERROR "art_add_dictionary: selection file '${selection_file}' does not exist")
#  endif()
#
#  # - Args may also be relative or absolute
#  foreach(_header ${ARGN})
#    if(NOT IS_ABSOLUTE "${_header}")
#      set(_header "${CMAKE_CURRENT_SOURCE_DIR}/${_header}")
#    endif()
#
#    if(NOT EXISTS "${_header}")
#      message(FATAL_ERROR "art_add_dictionary: input file '${_header}' does not exist")
#    endif()
#
#    list(APPEND headers "${_header}")
#  endforeach()
#
#  # Flags, include directories and definitions
#  art_get_genreflex_flags(__art_genreflex_flags)
#
#  # Any directory scope inc dirs
#  set(__art_genreflex_includes "-I${CMAKE_CURRENT_SOURCE_DIR}")
#  get_directory_property(local_incdirs INCLUDE_DIRECTORIES)
#  foreach(_incdir ${local_incdirs})
#    list(APPEND __art_genreflex_includes "-I${_incdir}")
#    message(STATUS "incdirs: ${_incdir}")
#  endforeach()
#  set(__art_genreflex_includes "-I$<JOIN:$<TARGET_PROPERTY:${__name},INCLUDE_DIRECTORIES>,\t-I>")
#
#  add_custom_command(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${__name}_dict.cpp
#           ${CMAKE_CURRENT_BINARY_DIR}/${__name}_map.cpp
#    COMMAND ${ROOT_genreflex_CMD}
#            ${headers}
#            -s ${selection_file}
#            #"-I$<JOIN:$<TARGET_PROPERTY:${__name},INCLUDE_DIRECTORIES>,\t-I>"
#            ${__art_genreflex_includes}
#            --capabilities=${__name}_map.cpp
#            -o ${__name}_dict.cpp
#    IMPLICIT_DEPENDS CXX ${headers}
#    DEPENDS ${selectionfile}
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#
#    )
#
#  # Equally, could follow Qt4 code generation style and just output
#  # paths to output sources...
#  add_library(${__name}_dict SHARED ${CMAKE_CURRENT_BINARY_DIR}/${__name}_dict.cpp)
#  add_library(${__name}_map  SHARED ${CMAKE_CURRENT_BINARY_DIR}/${__name}_map.cpp)
#endfunction()




