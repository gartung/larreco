art_add_module(BezierTrackerModule_module BezierTrackerModule_module.cc)

art_add_module(CCTrackMaker_module CCTrackMaker_module.cc)

art_add_module(CosmicTracker_module CosmicTracker_module.cc)

art_add_module(DumpTracks_module DumpTracks_module.cc)

art_add_module(FeatureTracker_module FeatureTracker_module.cc)

art_add_module(MagDriftAna_module MagDriftAna_module.cc)

art_add_module(SeedAna_module SeedAna_module.cc)

art_add_module(SeedFinderModule_module SeedFinderModule_module.cc)

art_add_module(SpacePointAna_module SpacePointAna_module.cc)

art_add_module(SpacePointCheater_module SpacePointCheater_module.cc)

art_add_module(SpacePointFinder_module SpacePointFinder_module.cc)

art_add_module(SpacePts_module SpacePts_module.cc)

art_add_module(Track3DKalmanHit_module Track3DKalmanHit_module.cc)

art_add_module(Track3DKalmanSPS_module Track3DKalmanSPS_module.cc)

art_add_module(Track3DKalman_module Track3DKalman_module.cc)

art_add_module(Track3Dreco_module Track3Dreco_module.cc)

art_add_module(TrackAna_module TrackAna_module.cc)

art_add_module(TrackCheater_module TrackCheater_module.cc)

art_add_module(TrackKalmanCheater_module TrackKalmanCheater_module.cc)

art_add_module(TrackStitcher_module TrackStitcher_module.cc)

install(TARGETS
     BezierTrackerModule_module
     CCTrackMaker_module
     CosmicTracker_module
     DumpTracks_module
     FeatureTracker_module
     MagDriftAna_module
     SeedAna_module
     SeedFinderModule_module
     SpacePointAna_module
     SpacePointCheater_module
     SpacePointFinder_module
     SpacePts_module
     Track3DKalmanHit_module
     Track3DKalmanSPS_module
     Track3DKalman_module
     Track3Dreco_module
     TrackAna_module
     TrackCheater_module
     TrackKalmanCheater_module
     TrackStitcher_module
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )


file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)
