IF(ALT_CMAKE)
INCLUDE(altCMakeLists.cmake)
ELSE()
art_make( LIBRARY_NAME  ClusterRecoUtil
          LIB_LIBRARIES
                        Utilities
                        ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
ENDIF()
