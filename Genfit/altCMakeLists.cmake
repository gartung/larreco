set(Genfit_HEADERS
     GFAbsBField.h
     GFAbsEnergyLoss.h
     GFAbsFinitePlane.h
     GFAbsGeoMatManager.h
     GFAbsRecoHit.h
     GFAbsTrackRep.h
     GFBookkeeping.h
     GFConstField.h
     GFDaf.h
     GFDetPlane.h
     GFEnergyLossBetheBloch.h
     GFEnergyLossBrems.h
     GFEnergyLossCoulomb.h
     GFException.h
     GFFieldManager.h
     GFGeoMatManager.h
     GFKalman.h
     GFMaterialEffects.h
     GFPlanarHitPolicy.h
     GFRecoHitFactory.h
     GFRecoHitIfc.h
     GFRecoHitProducer.h
     GFRectFinitePlane.h
     GFSpacepointHitPolicy.h
     GFTrack.h
     GFTrackCand.h
     GFWireHitPolicy.h
     GFWirepointHitPolicy.h
     GeaneMCApplication.h
     GeaneTrackRep2.h
     PointHit.h
     RKTrackRep.h
     SlTrackRep.h
     )

add_library(Genfit SHARED
     ${Genfit_HEADERS}
     GFAbsEnergyLoss.cxx
     GFAbsFinitePlane.cxx
     GFAbsGeoMatManager.cxx
     GFAbsRecoHit.cxx
     GFAbsTrackRep.cxx
     GFBookkeeping.cxx
     GFConstField.cxx
     GFDaf.cxx
     GFDetPlane.cxx
     GFEnergyLossBetheBloch.cxx
     GFEnergyLossBrems.cxx
     GFEnergyLossCoulomb.cxx
     GFException.cxx
     GFFieldManager.cxx
     GFGeoMatManager.cxx
     GFKalman.cxx
     GFMaterialEffects.cxx
     GFPlanarHitPolicy.cxx
     GFRecoHitFactory.cxx
     GFRecoHitProducer.cxx
     GFRectFinitePlane.cxx
     GFSpacepointHitPolicy.cxx
     GFTrack.cxx
     GFTrackCand.cxx
     GFWireHitPolicy.cxx
     GFWirepointHitPolicy.cxx
     GeaneMCApplication.cxx
     GeaneTrackRep2.cxx
     PointHit.cxx
     RKTrackRep.cxx
     SlTrackRep.cxx
     )

target_link_libraries(Genfit
     art::art_Framework_Core
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     art::art_Framework_Services_Registry
     FNALCore::FNALCore
     ${ROOT_BASIC_LIB_LIST}
     ${ROOT_GEOM}
     ${ROOT_EG}
     ${ROOT_VMC}
)

install(TARGETS
     Genfit
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime
     )

install(FILES ${Genfit_HEADERS} DESTINATION 
     ${CMAKE_INSTALL_INCLUDEDIR}/Genfit COMPONENT Development)

