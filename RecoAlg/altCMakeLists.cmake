set(HEADERS
     APAGeometryAlg.h
     BezierTrackerAlgorithm.h
     CCHitFinderAlg.h
     CCHitRefinerAlg.h
     ClusterCrawlerAlg.h
     ClusterMatchAlg.h
     ClusterMergeAlg.h
     ClusterMergeHelper.h
     ClusterParamsImportWrapper.h
     CornerFinderAlg.h
     DBScanAlg.h
     DisambigAlg.h
     EndPointAlg.h
     HoughBaseAlg.h
     KalmanFilterAlg.h
     LinFitAlg.h
     RootMathFunctor.h
     SeedFinderAlgorithm.h
     SmallClusterFinderAlg.h
     SpacePointAlg.h
     SpacePointAlg_TimeSort.h
     StitchAlg.h
     TrackLineFitAlg.h
     TrackMomentumCalculator.h
     TrackTrajectoryAlg.h
     fuzzyClusterAlg.h
     )

add_library(RecoAlg SHARED
     ${HEADERS}
     APAGeometryAlg.cxx
     BezierTrackerAlgorithm.cxx
     CCHitFinderAlg.cxx
     CCHitRefinerAlg.cxx
     ClusterCrawlerAlg.cxx
     ClusterMatchAlg.cxx
     ClusterMergeHelper.cxx
     CornerFinderAlg.cxx
     DBScanAlg.cxx
     DisambigAlg.cxx
     EndPointAlg.cxx
     HoughBaseAlg.cxx
     KalmanFilterAlg.cxx
     LinFitAlg.cxx
     SeedFinderAlgorithm.cxx
     SmallClusterFinderAlg.cxx
     SpacePointAlg.cxx
     SpacePointAlg_TimeSort.cxx
     StitchAlg.cxx
     TrackLineFitAlg.cxx
     TrackMomentumCalculator.cxx
     TrackTrajectoryAlg.cxx
     fuzzyClusterAlg.cxx
     )

target_link_libraries(RecoAlg
     larsoft::RecoBase
     larsoft::Simulation
     larsoft::Filters
     larsoft::Geometry
     larsoft::RecoObjects
     larsoft::AnalysisBase
     larsoft::Utilities
     art::art_Framework_Core
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     art::art_Framework_Services_Registry
     FNALCore::FNALCore
     ${ROOT_BASIC_LIB_LIST}
     ${ROOT_MINUIT}
     ${ROOT_MINUIT2}	
)

install(TARGETS
     RecoAlg
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

install(FILES ${HEADERS} DESTINATION 
     ${CMAKE_INSTALL_INCLUDEDIR}/RecoAlg COMPONENT Development)


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)

add_subdirectory(CMTool)
add_subdirectory(ClusterRecoUtil)
add_subdirectory(Cluster3DAlgs)
