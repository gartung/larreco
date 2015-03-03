set( PACKAGE CMTToolBase )


set( HEADERS
     CBoolAlgoBase.h
     CFloatAlgoBase.h
     CMAlgoBase.h
     CMManagerBase.h
     CMTException.h
     CMTool-TypeDef.h
     CMatchBookKeeper.h
     CMatchManager.h
     CMergeBookKeeper.h
     CMergeManager.h
     CPriorityAlgoBase.h
     )

set( SOURCES
     ${HEADERS}
     CBoolAlgoBase.cxx
     CFloatAlgoBase.cxx
     CMAlgoBase.cxx
     CMManagerBase.cxx
     CMTException.cxx
     CMatchBookKeeper.cxx
     CMatchManager.cxx
     CMergeBookKeeper.cxx
     CMergeManager.cxx
     CPriorityAlgoBase.cxx
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
     ${CMAKE_INSTALL_INCLUDEDIR}/RecoAlg/CMTool/CMTToolBase COMPONENT Development)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


