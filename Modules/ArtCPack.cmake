# - Cpackaging file for Art

#-----------------------------------------------------------------------
# Generic settings
#
# - Package name is project name...
# - Versioning
set(CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH})

#-----------------------------------------------------------------------
# Specifics for Source Package
#
set(CPACK_SOURCE_GENERATOR "TBZ2;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
  "${PROJECT_BINARY_DIR}"
  "/\\\\.git"
  "\\\\.swp$"
  )

#-----------------------------------------------------------------------
# Must allways include Cpack module **last**
#
include(CPack)

