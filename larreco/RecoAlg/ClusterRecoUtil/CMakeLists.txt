art_make( LIB_LIBRARIES
                        lardata_Utilities
                        larcorealg_Geometry
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_LIBRARIES}
        )

install_headers()
install_fhicl()
install_source()
