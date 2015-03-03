set( PACKAGE CMTAlgPriority )


set( HEADERS
     CMTAlgPriority-TypeDef.h
     CPAlgoArray.h
     CPAlgoIgnoreTracks.h
     CPAlgoNHits.h
     CPAlgoPolyArea.h
     CPAlgoQSum.h
     )

set( SOURCES
     ${HEADERS}
     CPAlgoArray.cxx
     CPAlgoIgnoreTracks.cxx
     CPAlgoNHits.cxx
     CPAlgoPolyArea.cxx
     CPAlgoQSum.cxx
     )

add_library(${PACKAGE} SHARED ${SOURCES})



target_link_libraries(${PACKAGE}
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
     ${PACKAGE}
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

install(FILES ${HEADERS} DESTINATION 
     ${CMAKE_INSTALL_INCLUDEDIR}/RecoAlg/CMTool/CMTAlgPriority COMPONENT Development)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


