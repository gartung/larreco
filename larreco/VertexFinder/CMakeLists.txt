set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation -Wno-maybe-uninitialized -Wno-overloaded-virtual")

art_make(
          EXCLUDE FeatureVertexFinder_module.cc CornerFinder_module.cc
          MODULE_LIBRARIES
                        larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        larevt_Filters
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        nusimdata_SimulationBase
                        lardata_RecoObjects
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_LIBRARIES}
                        ${ROOT_GenVector_LIBRARY}
                        larreco_RecoAlg
         )


simple_plugin ( FeatureVertexFinder "module"
                        larreco_RecoAlg
                        larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        larevt_Filters
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_LIBRARIES}
              )

simple_plugin ( CornerFinder "module"
                        larreco_RecoAlg
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_LIBRARIES}
              )

install_headers()
install_fhicl()
install_source()


