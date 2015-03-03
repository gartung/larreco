# macro for building ROOT dictionaries
#
# USAGE:
# build_dictionary( [<dictionary_name>]
#                   [COMPILE_FLAGS <flags>]
#                   [DICT_NAME_VAR <var>]
#                   [DICTIONARY_LIBRARIES <library list>]
#                   [NO_INSTALL]
#                   [DICT_FUNCTIONS])
#
# * <dictionary_name> defaults to a name based on the current source
# code subdirectory.
#
# * ${REFLEX} is always appended to the library list (even if it is
# empty).
#
# * Specify NO_INSTALL when building a dictionary for tests.
#
# * The default behavior is to generate a dictionary for data only. Use
# the DICT_FUNCTIONS option to reactivate the generation of dictionary
# entries for functions.
#
# * If DICT_NAME_VAR is specified, <var> will be set to contain the
# dictionary name.
#
# * Any other macros or functions in this file are for internal use
# only.
#
########################################################################
include(CMakeParseArguments)

# define flags for genreflex
set(GENREFLEX_FLAGS
  --deep
  --iocomments
  --fail_on_warnings
  --capabilities=classes_ids.cc
  --gccxmlopt=--gccxml-compiler
  --gccxmlopt=${CMAKE_CXX_COMPILER}
  -D_REENTRANT
  -DGNU_SOURCE
  -DGNU_GCC
  -DPROJECT_NAME="${PROJECT_NAME}"
  -DPROJECT_VERSION="${${PROJECT_NAME}_VERSION}"
  -D__STRICT_ANSI__
  )

#-----------------------------------------------------------------------
# macro _set_dictionary_name()
#
macro(_set_dictionary_name)
  if(PACKAGE_TOP_DIRECTORY)
    string(REGEX REPLACE "^${PACKAGE_TOP_DIRECTORY}/(.*)" "\\1" CURRENT_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}")
  else()
    string(REGEX REPLACE "^${CMAKE_SOURCE_DIR}/(.*)" "\\1" CURRENT_SUBDIR "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  string(REGEX REPLACE "/" "_" dictname "${CURRENT_SUBDIR}" )
endmacro()

#-----------------------------------------------------------------------
# function _generate_dictionary()
#
function(_generate_dictionary)
  cmake_parse_arguments(GD
    "DICT_FUNCTIONS"
    "NAME"
    "DICTIONARY_LIBRARIES"
    ${ARGN}
    )
  set(generate_dictionary_usage "_generate_dictionary( [DICT_FUNCTIONS] [DICTIONARY_LIBRARIES lib1 lib2 ... NAME [dictionary_name] )")

  # Error on unknown arguments
  if(NOT ${GD_NAME})
    _set_dictionary_name()
  else()
    set(dictname ${GD_NAME})
    list(LENGTH GD_UNPARSED_ARGUMENTS n_bad_args)
    if(n_bad_args GREATER 1)
      list(REMOVE_AT GD_UNPARSED_ARGUMENTS 0)
      message("_GENERATE_DICTIONARY: unwanted extra arguments: ${GD_UNPARSED_ARGUMENTS}")
      message(SEND_ERROR ${generate_dictionary_usage})
    endif()
  endif()

  # Configure genreflex command line args
  if(NOT GD_DICT_FUNCTIONS AND NOT CET_DICT_FUNCTIONS)
    set(GENREFLEX_FLAGS ${GENREFLEX_FLAGS} --dataonly)
  endif()

  # Local Include dirs
  get_directory_property(genpath INCLUDE_DIRECTORIES)
  foreach(inc ${genpath})
    set(GENREFLEX_INCLUDES ${GENREFLEX_INCLUDES} -I ${inc})
  endforeach()

  # Local compile defs
  get_directory_property(compile_defs COMPILE_DEFINITIONS)
  foreach(def ${compile_defs})
    set(GENREFLEX_FLAGS ${GENREFLEX_FLAGS} -D${def})
  endforeach()

  # Include dirs and compile defs specific to requested link libs
  foreach(extralib ${GD_DICTIONARY_LIBRARIES})
    if(TARGET ${extralib})
      # Yeah, these are knarly - note the following:
      # - They are well formed, and taken from
      # http://www.cmake.org/cmake/help/v3.0/manual/cmake-generator-expressions.7.html
      # - The heavy wrapping is needed because we don't want to
      #   pass a blank '-I' or '-D' to genreflex - it really doesn't
      #   like empty -Ds
      # - A '\t' is used to join lists to prevent escaping when
      #  the expression is expanded in the custom_command
      #
      # - Include dirs
      set(GENREFLEX_INCLUDES ${GENREFLEX_INCLUDES} "$<$<BOOL:$<TARGET_PROPERTY:${extralib},INCLUDE_DIRECTORIES>>:-I$<JOIN:$<TARGET_PROPERTY:${extralib},INCLUDE_DIRECTORIES>,\t-I>>")

      # - Compile defs
      set(GENREFLEX_FLAGS ${GENREFLEX_FLAGS} "$<$<BOOL:$<TARGET_PROPERTY:${extralib},COMPILE_DEFINITIONS>>:-D$<JOIN:$<TARGET_PROPERTY:${extralib},COMPILE_DEFINITIONS>,\t-D>>")
    endif()
  endforeach()

  if(${GENREFLEX_CLEANUP} MATCHES "TRUE")
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_dict.cpp
             ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_map.cpp
      COMMAND ${ROOT_genreflex_CMD}
              ${CMAKE_CURRENT_SOURCE_DIR}/classes.h
          	  -s ${CMAKE_CURRENT_SOURCE_DIR}/classes_def.xml
		          -I ${CMAKE_SOURCE_DIR}
		          -I ${CMAKE_CURRENT_SOURCE_DIR}
		          ${GENREFLEX_INCLUDES}
              ${GENREFLEX_FLAGS}
        	    -o ${dictname}_dict.cpp || {rm -f ${dictname}_dict.cpp\; /bin/false\;}
      COMMAND ${CMAKE_COMMAND} -E copy classes_ids.cc ${dictname}_map.cpp
      COMMAND ${CMAKE_COMMAND} -E remove -f classes_ids.cc
      IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/classes.h
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/classes_def.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
  else()
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_dict.cpp
             ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_map.cpp
      COMMAND ${ROOT_genreflex_CMD}
              ${CMAKE_CURRENT_SOURCE_DIR}/classes.h
           	  -s ${CMAKE_CURRENT_SOURCE_DIR}/classes_def.xml
		          -I ${CMAKE_SOURCE_DIR}
		          -I ${CMAKE_CURRENT_SOURCE_DIR}
		          ${GENREFLEX_INCLUDES}
              ${GENREFLEX_FLAGS}
        	    -o ${dictname}_dict.cpp
      COMMAND ${CMAKE_COMMAND} -E copy classes_ids.cc ${dictname}_map.cpp
      COMMAND ${CMAKE_COMMAND} -E remove -f classes_ids.cc
      IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/classes.h
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/classes_def.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
  endif()

  set_source_files_properties(
    ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_dict.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_map.cpp
    PROPERTIES GENERATED 1
    )
endfunction()

#-----------------------------------------------------------------------
# function build_dictionary()
#
function(build_dictionary)
  set(build_dictionary_usage "USAGE: build_dictionary( [dictionary_name] [DICTIONARY_LIBRARIES <library list>] [COMPILE_FLAGS <flags>] [DICT_NAME_VAR <var>] [NO_INSTALL] )")

  cmake_parse_arguments(BD
    "NOINSTALL;NO_INSTALL;DICT_FUNCTIONS"
    "DICT_NAME_VAR"
    "DICTIONARY_LIBRARIES;COMPILE_FLAGS"
    ${ARGN}
    )

  # Handle obsolete and unknown args
  if(BD_NOINSTALL OR BD_NO_INSTALL)
    message(SEND_ERROR "build_dictionary no longer performs any install operation")
  endif()

  if(BD_UNPARSED_ARGUMENTS)
    list(LENGTH BD_UNPARSED_ARGUMENTS dlen)
    if(dlen GREATER 1 )
	    message("BUILD_DICTIONARY: Too many arguments. ${ARGV}")
	    message(SEND_ERROR  ${build_dictionary_usage})
    endif()
    list(GET BD_UNPARSED_ARGUMENTS 0 dictname)
  else()
    _set_dictionary_name()
  endif()

  if(BD_DICT_NAME_VAR)
    set(${BD_DICT_NAME_VAR} ${dictname} PARENT_SCOPE)
  endif()

  if(BD_DICTIONARY_LIBRARIES)
    set(dictionary_liblist ${BD_DICTIONARY_LIBRARIES})
  endif()

  list(APPEND dictionary_liblist ${ROOT_Core_LIBRARY} ${ROOT_Reflex_LIBRARY})
  _generate_dictionary(DICTIONARY_LIBRARIES ${dictionary_liblist} NAME ${dictname})

  add_library(${dictname}_dict SHARED ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_dict.cpp )
  add_library(${dictname}_map  SHARED ${CMAKE_CURRENT_BINARY_DIR}/${dictname}_map.cpp )

  if(BD_COMPILE_FLAGS)
    set_target_properties(${dictname}_dict ${dictname}_map
      PROPERTIES COMPILE_FLAGS ${BD_COMPILE_FLAGS})
  endif()
  target_link_libraries(${dictname}_dict ${dictionary_liblist})
  target_link_libraries(${dictname}_map  ${dictionary_liblist})
  add_dependencies(${dictname}_map ${dictname}_dict)
endfunction()
