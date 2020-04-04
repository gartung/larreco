include(${Geant4_USE_FILE})
include_directories ( ${Geant4_INCLUDE_DIRS} )

art_make(MODULE_LIBRARIES
           larreco_RecoAlg
           ${MF_MESSAGELOGGER}
         )
install_headers()
install_fhicl()
