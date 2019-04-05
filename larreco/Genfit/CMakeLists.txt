
art_make( LIB_LIBRARIES
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_LIBRARIES}
                        ${ROOT_Geom_LIBRARY}
                        ${ROOT_EG_LIBRARY}
                        ${ROOT_VMC_LIBRARY}
      )

install_headers()
install_fhicl()
install_source()
