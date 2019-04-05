include(${Geant4_USE_FILE})
include_directories ( ${Geant4_INCLUDE_DIRS} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS} -Wno-shadow -Wno-unused-parameter")

art_make(
          MODULE_LIBRARIES  larreco_RecoAlg
             lardata_RecoObjects larreco_Deprecated
             lardataobj_RecoBase
             larcorealg_Geometry
             larcore_Geometry_Geometry_service
             lardata_Utilities
             nutools_MagneticField_MagneticField_service
             nusimdata_SimulationBase
             ${ART_FRAMEWORK_CORE}
             ${ART_FRAMEWORK_PRINCIPAL}
             ${ART_FRAMEWORK_SERVICES_REGISTRY}
             ${ART_FRAMEWORK_SERVICES_BASIC}
             ${ART_ROOT_IO_TFILE_SUPPORT}
             ${ART_ROOT_IO_TFILESERVICE_SERVICE}
             ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
             art_Persistency_Common
             art_Utilities
             canvas
             ${MF_MESSAGELOGGER}
             ${FHICLCPP}
             cetlib cetlib_except
             ${CLHEP_LIBRARIES}
             ${ROOT_LIBRARIES}
             ${Geant4_LIBRARIES}
         )
install_headers()
install_fhicl()
install_source()


