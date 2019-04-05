include(${Geant4_USE_FILE})
include_directories ( ${Geant4_INCLUDE_DIRS} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Geant4_CXX_FLAGS} -Wno-shadow -Wno-unused-parameter -Wno-ignored-qualifiers")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation -Wno-maybe-uninitialized -Wno-overloaded-virtual")

art_make(
          EXCLUDE Track3DKalman_module.cc Track3DKalmanSPS_module.cc DumpTracks_module.cc KalmanFilterFitTrackMaker_tool.cc
          MODULE_LIBRARIES  larreco_RecoAlg
             larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
             lardata_RecoObjects larreco_Deprecated
             lardataobj_RecoBase
             larcorealg_Geometry
             larcore_Geometry_Geometry_service
             larreco_Calorimetry
             lardataobj_AnalysisBase
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

simple_plugin( Track3DKalman "module"
                 larreco_Genfit
                 larreco_RecoAlg
                 larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                 lardata_RecoObjects larreco_Deprecated
                 lardataobj_RecoBase
                 larcorealg_Geometry
                 larcore_Geometry_Geometry_service
                 larreco_Calorimetry
                 lardataobj_AnalysisBase
                 lardata_Utilities
                 nutools_RandomUtils_NuRandomService_service
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

        )

simple_plugin( Track3DKalmanSPS "module"
          larreco_Genfit
          larreco_RecoAlg
          larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
          lardata_RecoObjects larreco_Deprecated
          lardataobj_RecoBase
          larcorealg_Geometry
          larcore_Geometry_Geometry_service
          larreco_Calorimetry
          lardataobj_AnalysisBase
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

       )

simple_plugin( KalmanFilterFitTrackMaker "tool"
             ${Boost_SYSTEM_LIBRARY}
             ${Boost_FILESYSTEM_LIBRARY}
             ${ROOT_Core_LIBRARY}
             canvas
             lardataobj_RecoBase
             larreco_RecoAlg
       )

install_headers()
install_fhicl()
install_source()


