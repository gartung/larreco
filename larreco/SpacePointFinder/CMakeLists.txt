include(${Geant4_USE_FILE})
include_directories ( ${Geant4_INCLUDE_DIRS} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS} -Wno-shadow -Wno-unused-parameter")

art_make(MODULE_LIBRARIES
           larreco_RecoAlg
           ${MF_MESSAGELOGGER}
        )
install_headers()
install_fhicl()
install_source()
