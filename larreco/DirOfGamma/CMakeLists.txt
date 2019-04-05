
art_make( 
          LIB_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larreco_RecoAlg
                        larreco_RecoAlg_PMAlg
                        lardata_Utilities
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_LIBRARIES}
                        ${Boost_SYSTEM_LIBRARY}
     MODULE_LIBRARIES larreco_RecoAlg
                      larreco_RecoAlg_PMAlg
                      larreco_DirOfGamma
                      larsim_MCCheater_BackTrackerService_service
                      larsim_MCCheater_ParticleInventoryService_service
                      lardataobj_RecoBase
                      larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
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
                      art_Persistency_Provenance
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
