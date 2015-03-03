
set(HEADERS
     DBScanAlg.h
     HoughSeedFinderAlg.h
     PCASeedFinderAlg.h
     ParallelHitsSeedFinderAlg.h
     PrincipalComponentsAlg.h
     SeedFinderAlgBase.h
     SkeletonAlg.h
     )

add_library(Cluster3DAlgs SHARED
     ${HEADERS}
     DBScanAlg.cxx
     HoughSeedFinderAlg.cxx
     PCASeedFinderAlg.cxx
     ParallelHitsSeedFinderAlg.cxx
     PrincipalComponentsAlg.cxx
     SkeletonAlg.cxx
     )

target_link_libraries(Cluster3DAlgs
     larsoft::Utilities
     art::art_Framework_Core
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     art::art_Framework_Services_Registry
     FNALCore::FNALCore
     ${ROOT_BASIC_LIB_LIST}
)

install(TARGETS
     Cluster3DAlgs
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

install(FILES ${HEADERS} DESTINATION 
     ${CMAKE_INSTALL_INCLUDEDIR}/RecoAlg/Cluster3DAlgs COMPONENT Development)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


