# ======================================================================
#  larreco main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larreco/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 3.10 or later
cmake_minimum_required (VERSION 3.10)

project(larreco VERSION 08.00.00 LANGUAGES CXX)

set(fw_dir config_data)
set(fcl_dir job)
set(gdml_dir gdml)
find_package( cetmodules )
include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_package( ROOT REQUIRED COMPONENTS GenVector Geom Minuit Minuit2 XMLIO Gdml TMVA EG VMC)
find_package( larsim REQUIRED  )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# not ready for c2 yet
if( DEFINED TENSORFLOW_DIR )
  # find tensorflow library
  find_package( Protobuf )
  find_package( Tensorflow )
  cet_find_library(TENSORFLOW NAMES tensorflow-core )
  cet_find_library(PROTOBUF NAMES protobuf )
  # This is a hack for a broken tensorflow link, which
  # requires -lrt for glibc earlier than 2.17 (i.e. SLF6)
  # We take advantage of OSTYPE defined by cetbuildtools
  if ( "${OSTYPE}" MATCHES "slf6" )
    set(FOR_TENSORFLOW "-lrt")
  endif ()
endif ()

# source
add_subdirectory(larreco)

# tests
add_subdirectory(test)
add_subdirectory(cmakeConfig)
# ups - table and config files
#add_subdirectory(ups)

# packaging utility
#include(UseCPack)
