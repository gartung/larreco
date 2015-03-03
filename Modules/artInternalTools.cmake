# - Developer tools for use in build of Art only

#-----------------------------------------------------------------------
# Copyright 2014 Ben Morgan <Ben.Morgan@warwick.ac.uk>
# Copyright 2014 University of Warwick

# - Needed core CMake support
include(GNUInstallDirs)

# - Standard arguments to pass to install(TARGETS ...) command
set(art_TARGET_INSTALL_ARGS
  EXPORT ArtLibraries
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT Runtime
  )

#-----------------------------------------------------------------------
# function art_set_standard_target_properties(<tgt>)
#          Apply standard global target properties to target <tgt>
#          Motivation - there are always some properties to apply
#                       globally, e.g. soversions
function(art_set_standard_target_properties __art_target)
  if(NOT TARGET ${__art_target})
    message(FATAL_ERROR "${__art_target} is not a target")
  endif()

  # If we can't see the version numbers, barf
  if((NOT art_VERSION) OR (NOT art_SOVERSION))
    message(FATAL_ERROR "Cannot access art's VERSION/SOVERSION")
  endif()


  set_target_properties(${__art_target}
    PROPERTIES
     VERSION ${art_VERSION}
     SOVERSION ${art_SOVERSION}
    )
endfunction()

