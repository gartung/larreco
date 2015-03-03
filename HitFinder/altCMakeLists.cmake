add_library(HitFinder SHARED
     HitAnaAlg.cxx
     HitAnaAlg.h
)

target_link_libraries(HitFinder
     larsoft::RecoBase
     art::art_Framework_Core
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     art::art_Framework_Services_Registry
     FNALCore::FNALCore
     )

art_add_module(APAHitFinder_module APAHitFinder_module.cc)

art_add_module(DisambigCheater_module DisambigCheater_module.cc)

art_add_module(FFTHitFinder_module FFTHitFinder_module.cc)

art_add_module(GausHitFinderAna_module GausHitFinderAna_module.cc)

art_add_module(GausHitFinder_module GausHitFinder_module.cc)

art_add_module(HitAnaModule_module HitAnaModule_module.cc)

art_add_module(HitCheater_module HitCheater_module.cc)

art_add_module(HitFinderAna_module HitFinderAna_module.cc)

art_add_module(MCHitAnaExample_module MCHitAnaExample_module.cc)

art_add_module(MCHitFinder_module MCHitFinder_module.cc)

art_add_module(RFFHitFinder_module RFFHitFinder_module.cc)

art_add_module(TTHitFinder_module TTHitFinder_module.cc)

install(TARGETS
     HitFinder
     APAHitFinder_module
     DisambigCheater_module
     FFTHitFinder_module
     GausHitFinderAna_module
     GausHitFinder_module
     HitAnaModule_module
     HitCheater_module
     HitFinderAna_module
     MCHitAnaExample_module
     MCHitFinder_module
     RFFHitFinder_module
     TTHitFinder_module
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )


install(FILES HitAnaAlg.h DESTINATION
     ${CMAKE_INSTALL_INCLUDEDIR}/HitFinder COMPONENT Development)

file(GLOB FHICL_FILES 
     [^.]*.fcl
)

install(FILES ${FHICL_FILES} DESTINATION job COMPONENT Runtime)

