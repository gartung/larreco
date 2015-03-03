
set( PACKAGE CMTAlgMatch )


set( HEADERS
     CFAlgo3DAngle.h
     CFAlgoArray.h
     CFAlgoChargeDistrib.h
     CFAlgoQRatio.h
     CFAlgoStartPointCompat.h
     CFAlgoStartPointMatch.h
     CFAlgoStartTimeCompat.h
     CFAlgoTimeOverlap.h
     CFAlgoTimeProf.h
     CFAlgoVolumeOverlap.h
     CFAlgoWireOverlap.h
     CFAlgoZOverlap.h
     CMTAlgMatch-TypeDef.h
     )

set( SOURCES
     ${HEADERS}
     CFAlgo3DAngle.cxx
     CFAlgoArray.cxx
     CFAlgoChargeDistrib.cxx
     CFAlgoQRatio.cxx
     CFAlgoStartPointCompat.cxx
     CFAlgoStartPointMatch.cxx
     CFAlgoStartTimeCompat.cxx
     CFAlgoTimeOverlap.cxx
     CFAlgoTimeProf.cxx
     CFAlgoVolumeOverlap.cxx
     CFAlgoWireOverlap.cxx
     CFAlgoZOverlap.cxx
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
     ${CMAKE_INSTALL_INCLUDEDIR}/RecoAlg/CMTool/CMTAlgMatch COMPONENT Development)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


