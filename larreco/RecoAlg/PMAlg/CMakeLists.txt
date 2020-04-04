art_make(LIB_LIBRARIES
           larcorealg_Geometry
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ROOT::Core
           ROOT::Matrix
           ROOT::Physics
           ROOT::Minuit2
           canvas
           ${MF_MESSAGELOGGER}
           cetlib_except
        )

install_headers()
install_fhicl()
