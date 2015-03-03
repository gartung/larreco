#string(REPLACE "/" "+" _ccv_target "${EXECUTABLE_OUTPUT_PATH}/checkClassVersion")

# - Use a property to configure the dynamic loader path for checkClassVersion
# It's done this way because we don't want to hard code anything in the
# actual macro or provide an argument (because it's used internally by ArtDictionary
# and that can be used by clients...
function(checkclassversion_append_path _path)
  set_property(GLOBAL
    APPEND
    PROPERTY CHECKCLASSVERSION_DYNAMIC_PATH ${_path}
    )
endfunction()


MACRO(check_class_version)
  CMAKE_PARSE_ARGUMENTS(ART_CCV
    "LIBRARIES"
    "UPDATE_IN_PLACE"
    ${ARGN}
    )
  IF(ART_CCV_LIBRARIES)
    MESSAGE(FATAL_ERROR "LIBRARIES option not supported at this time: "
      "ensure your library is linked to any necessary libraries not already pulled in by ART.")
  ENDIF()
  IF(ART_CCV_UPDATE_IN_PLACE)
    SET(ART_CCV_EXTRA_ARGS ${ART_CCV_EXTRA_ARGS} "-G")
  ENDIF()
  IF(NOT dictname)
    MESSAGE(FATAL_ERROR "CHECK_CLASS_VERSION must be called after BUILD_DICTIONARY.")
  ENDIF()
  IF(ROOT_python_FOUND)
    # Use dynamic path lookup if configured
    get_property(ART_CCV_DYNAMIC_PATH GLOBAL PROPERTY CHECKCLASSVERSION_DYNAMIC_PATH)
    set(ART_CCV_DYNAMIC_PATH_ARGS)
    foreach(_path ${ART_CCV_DYNAMIC_PATH})
      list(APPEND ART_CCV_DYNAMIC_PATH_ARGS -L${_path})
    endforeach()

    # Add the check to the end of the dictionary building step.
    add_custom_command(TARGET ${dictname}_dict POST_BUILD
      COMMAND art::checkClassVersion
      ${ART_CCV_EXTRA_ARGS}
      ${ART_CCV_DYNAMIC_PATH_ARGS}
      -l$<TARGET_FILE:${dictname}_dict>
      -x ${CMAKE_CURRENT_SOURCE_DIR}/classes_def.xml
      COMMENT "Running checkClassVersion on ${dictname}_dict"
      VERBATIM
      )
    if (NOT TARGET art::art_Framework_Core)
      # If we're in art, we need to be sure that CheckClassVersion and
      # art_Framework_Core are already built; if we're outside art, this
      # is a given.
      add_dependencies(${dictname}_dict art::checkClassVersion art_Framework_Core)
    endif()
  ENDIF()
ENDMACRO()
