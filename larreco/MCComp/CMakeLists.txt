#add_subdirectory(job)

art_make(
          LIB_LIBRARIES    larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                           lardataobj_RecoBase
                           lardata_Utilities
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
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
                           ${ROOT_Geom_LIBRARY}
                           ${ROOT_XMLIO_LIBRARY}
                           ${ROOT_Gdml_LIBRARY}
                           ${ROOT_LIBRARIES}
          MODULE_LIBRARIES larreco_MCComp
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_BASIC}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${FHICLCPP}
                           cetlib cetlib_except
                           ${ROOT_Geom_LIBRARY}
                           ${ROOT_XMLIO_LIBRARY}
                           ${ROOT_Gdml_LIBRARY}
                           ${ROOT_LIBRARIES}
         )

# install_headers()
install_fhicl()
install_source()
