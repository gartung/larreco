art_make( LIB_LIBRARIES larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        lardataobj_RecoBase
                        larreco_Calorimetry
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${CLHEP_LIBRARIES}
                        ${ROOT_LIBRARIES}
                        ${ROOT_Minuit_LIBRARY}
                        ${ROOT_Minuit2_LIBRARY}
                        ${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_fhicl()
install_source()
