# - Try to find SQLite3
# Once done this will define
#
#  SQLite3_FOUND - system has Sqlite
#  SQLite3_INCLUDE_DIRS - the Sqlite include directory
#  SQLite3_LIBRARIES - Link these to use Sqlite
#  SQLite3_DEFINITIONS - Compiler switches required for using Sqlite
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

# Copyright (c) 2014, Ben Morgan, <Ben.Morgan@warwick.ac.uk>
# Copyright (c) 2008, Gilles Caulier, <caulier.gilles@gmail.com>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if(SQLite3_INCLUDE_DIR AND SQLite3_LIBRARY)
   # in cache already
   SET(SQLite3_FIND_QUIETLY TRUE)
endif()

# Optionally use pkg-config to get the directories and then use these values
# in the FIND_PATH() and FIND_LIBRARY() calls
if(NOT WIN32)
  find_package(PkgConfig)
  pkg_check_modules(PC_SQLite3 QUIET sqlite3)
  set(SQLite3_DEFINITIONS ${PC_SQLite3_CFLAGS_OTHER})
endif()

# - Find header
find_path(SQLite3_INCLUDE_DIR NAMES sqlite3.h
  PATHS
  ${PC_SQLite3_INCLUDEDIR}
  ${PC_SQLite3_INCLUDE_DIRS}
)

# - Extract version
if(SQLite3_INCLUDE_DIR AND EXISTS "${SQLite3_INCLUDE_DIR}/sqlite3.h")
  file(STRINGS "${SQLite3_INCLUDE_DIR}/sqlite3.h" SQLite3_H REGEX "^#define SQLITE_VERSION *\"[^\"]*\"$")

  string(REGEX REPLACE "^.*SQLITE_VERSION *\"([0-9]+).*$" "\\1" SQLite3_VERSION_MAJOR "${SQLite3_H}")
  string(REGEX REPLACE "^.*SQLITE_VERSION *\"[0-9]+\\.([0-9]+).*$" "\\1" SQLite3_VERSION_MINOR  "${SQLite3_H}")
  string(REGEX REPLACE "^.*SQLITE_VERSION *\"[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" SQLite3_VERSION_PATCH "${SQLite3_H}")
  set(SQLite3_VERSION_STRING "${SQLite3_VERSION_MAJOR}.${SQLite3_VERSION_MINOR}.${SQLite3_VERSION_PATCH}")

  # only append a TWEAK version if it exists:
  set(SQLite3_VERSION_TWEAK "")
  if( "${SQLite3_H}" MATCHES "^.*SQLITE_VERSION *\"[0-9]+\\.[0-9]+\\.[0-9]+\\.([0-9]+).*$")
    set(SQLite3_VERSION_TWEAK "${CMAKE_MATCH_1}")
    set(SQLite3_VERSION_STRING "${SQLite3_VERSION_STRING}.${SQLite3_VERSION_TWEAK}")
  endif()

  set(SQLite3_MAJOR_VERSION "${SQLite3_VERSION_MAJOR}")
  set(SQLite3_MINOR_VERSION "${SQLite3_VERSION_MINOR}")
  set(SQLite3_PATCH_VERSION "${SQLite3_VERSION_PATCH}")
endif()


find_library(SQLite3_LIBRARY NAMES sqlite3
  PATHS
  ${PC_SQLite3_LIBDIR}
  ${PC_SQLite3_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SQLite3
  FOUND_VAR SQLite3_FOUND
  REQUIRED_VARS SQLite3_INCLUDE_DIR SQLite3_LIBRARY
  VERSION_VAR SQLite3_VERSION_STRING
  )

# show the SQLite_INCLUDE_DIR and SQLite_LIBRARIES variables only in the advanced view
mark_as_advanced(SQLite3_INCLUDE_DIR SQLite3_LIBRARY)

# Final setup
if(SQLite3_FOUND)
  set(SQLite3_INCLUDE_DIRS ${SQLite3_INCLUDE_DIR})
  set(SQLite3_LIBRARIES ${SQLite3_LIBRARY})
endif()

