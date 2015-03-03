########################################################################
# art_add_dictionary
#
# Wrapper around cetbuildtools' build_dictionary featuring the addition
# of commonly required libraries to the dictionary library link list,
# and the use of the check_class_version to update checksums and
# class versions for dictionary items.
#
####################################
# Options and Arguments
#
# UPDATE_IN_PLACE
#   Passed through to check_class_version.
#
# DICT_FUNCTIONS
#   Passed through to build_dictionary.
#
# DICT_NAME_VAR
#   Passed through to build_dictionary.
#
# DICTIONARY_LIBRARIES
#   Passed through to build_dictionary with additions.
#
# COMPILE_FLAGS
#   Passed through to build_dictionary.
#
# USE_PRODUCT_NAME
#   Passed through to build_dictionary.
#
#########################################################################
include(BuildDictionary)
include(CMakeParseArguments)
include(CheckClassVersion)

function(art_dictionary)
  message(WARNING "art_dictionary is deprecated, use art_add_dictionary")
  art_add_dictionary(${ARGN})
endfunction()

function(art_add_dictionary)
  cmake_parse_arguments(AD
    "UPDATE_IN_PLACE;DICT_FUNCTIONS;USE_PRODUCT_NAME"
    "DICT_NAME_VAR"
    "DICTIONARY_LIBRARIES;COMPILE_FLAGS"
    ${ARGN}
    )

  # Setup common libs required for linking
  # We can use target names consistently across build and client because
  # of import/export of targets
  set(AD_DICTIONARY_LIBRARIES
    ${art_IMPORT_NAMESPACE}art_Persistency_Common
    ${art_IMPORT_NAMESPACE}art_Utilities
    FNALCore::FNALCore
    ${AD_DICTIONARY_LIBRARIES}
    )

  build_dictionary(DICT_NAME_VAR dictname
    DICTIONARY_LIBRARIES ${AD_DICTIONARY_LIBRARIES}
    ${AD_UNPARSED_ARGUMENTS}
    ${extra_args})

  # "returns"
  # We *probably* don't care about this as it only
  # appears to be relevant when installing source
  # code, and we almost certainly don't want to do that
  # for generated code.
  if (cet_generated_code) # Bubble up to top scope.
    set(cet_generated_code ${cet_generated_code} PARENT_SCOPE)
  endif()

  if (AD_DICT_NAME_VAR)
    set (${AD_DICT_NAME_VAR} ${dictname} PARENT_SCOPE)
  endif()

  if(AD_UPDATE_IN_PLACE)
    set(AD_CCV_ARGS ${AD_CCV_ARGS} "UPDATE_IN_PLACE" ${AD_UPDATE_IN_PLACE})
  endif()

  check_class_version(${AD_LIBRARIES} UPDATE_IN_PLACE ${AD_CCV_ARGS})
endfunction()
