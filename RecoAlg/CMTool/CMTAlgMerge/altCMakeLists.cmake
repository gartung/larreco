
set( PACKAGE CMTAlgMerge )

set( HEADERS
     CBAlgoAngleAlign.h
     CBAlgoAngleCompat.h
     CBAlgoAngleIncompat.h
     CBAlgoAngleSeparate.h
     CBAlgoArray.h
     CBAlgoCenterOfMass.h
     CBAlgoCenterOfMassSmall.h
     CBAlgoFake.h
     CBAlgoMergeAll.h
     CBAlgoMergeTinyWithBig.h
     CBAlgoOutOfConeSeparate.h
     CBAlgoPolyContain.h
     CBAlgoPolyHitOverlap.h
     CBAlgoPolyOverlap.h
     CBAlgoPolyShortestDist.h
     CBAlgoProhibitAllTracks.h
     CBAlgoProhibitBigClusters.h
     CBAlgoShortestDist.h
     CBAlgoShortestDistNonEndPoint.h
     CBAlgoShortestDistSmallCluster.h
     CBAlgoStartInCone.h
     CBAlgoStartInPoly.h
     CBAlgoStartNearEnd.h
     CBAlgoStartTrack.h
     CBAlgoTrackSeparate.h
     CMTAlgMerge-TypeDef.h
     )

set( SOURCES
     ${HEADERS}
     CBAlgoAngleAlign.cxx
     CBAlgoAngleCompat.cxx
     CBAlgoAngleIncompat.cxx
     CBAlgoAngleSeparate.cxx
     CBAlgoArray.cxx
     CBAlgoCenterOfMass.cxx
     CBAlgoCenterOfMassSmall.cxx
     CBAlgoFake.cxx
     CBAlgoMergeAll.cxx
     CBAlgoMergeTinyWithBig.cxx
     CBAlgoOutOfConeSeparate.cxx
     CBAlgoPolyContain.cxx
     CBAlgoPolyHitOverlap.cxx
     CBAlgoPolyOverlap.cxx
     CBAlgoPolyShortestDist.cxx
     CBAlgoProhibitAllTracks.cxx
     CBAlgoProhibitBigClusters.cxx
     CBAlgoShortestDist.cxx
     CBAlgoShortestDistNonEndPoint.cxx
     CBAlgoShortestDistSmallCluster.cxx
     CBAlgoStartInCone.cxx
     CBAlgoStartInPoly.cxx
     CBAlgoStartNearEnd.cxx
     CBAlgoStartTrack.cxx
     CBAlgoTrackSeparate.cxx
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
     ${CMAKE_INSTALL_INCLUDEDIR}/RecoAlg/CMTool/CMTAlgMerge COMPONENT Development)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)


