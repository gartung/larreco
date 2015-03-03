
set(HEADERS
     CRUException.h
     ClusterParams.h
     ClusterParamsAlg.h
     ClusterParamsAlgBase.h
     LazyClusterParamsAlg.h
     OverriddenClusterParamsAlg.h
     Polygon2D.h
     StandardClusterParamsAlg.h
     )

add_library(ClusterRecoUtil SHARED
     ${HEADERS}
     CRUException.cxx
     ClusterParamsAlg.cxx
     LazyClusterParamsAlg.cxx
     Polygon2D.cxx
     StandardClusterParamsAlg.cxx
     )

target_link_libraries(ClusterRecoUtil
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
     ClusterRecoUtil
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

install(FILES ${HEADERS} DESTINATION 
     ${CMAKE_INSTALL_INCLUDEDIR}/RecoAlg/ClusterRecoUtil COMPONENT Development)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)

