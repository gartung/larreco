art_make( LIBRARY_NAME  ClusterRecoUtilAlg
          LIB_LIBRARIES
                        ClusterRecoUtil
                        Utilities
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
