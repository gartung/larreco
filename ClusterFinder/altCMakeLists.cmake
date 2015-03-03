add_library(ClusterFinder SHARED
     ClusterCreator.h
     ClusterCreator.cxx
)

target_link_libraries(ClusterFinder
     larsoft::Geometry
     art::art_Framework_Core
     art::art_Framework_Principal
     art::art_Persistency_Provenance
     art::art_Utilities
     art::art_Framework_Services_Registry
     FNALCore::FNALCore
)

install(FILES ClusterCreator.h DESTINATION 
     ${CMAKE_INSTALL_INCLUDEDIR}/ClusterFinder COMPONENT Development)

art_add_module(ClusterAna_module ClusterAna_module.cc)

art_add_module(ClusterCheater_module ClusterCheater_module.cc)

art_add_module(ClusterCrawler_module ClusterCrawler_module.cc)

art_add_module(ClusterPCA_module ClusterPCA_module.cc)

art_add_module(Cluster3D_module Cluster3D_module.cc)

art_add_module(DBclusterAna_module DBclusterAna_module.cc)

art_add_module(DBcluster_module DBcluster_module.cc)

art_add_module(DumpClusters_module DumpClusters_module.cc)

art_add_module(EndPointModule_module EndPointModule_module.cc)

art_add_module(HoughLineFinderAna_module HoughLineFinderAna_module.cc)

art_add_module(HoughLineFinder_module HoughLineFinder_module.cc)

art_add_module(LineMerger_module LineMerger_module.cc)

art_add_module(SimpleClusterMerger_module SimpleClusterMerger_module.cc)

art_add_module(SmallClusterFinder_module SmallClusterFinder_module.cc)

art_add_module(fuzzyCluster_module fuzzyCluster_module.cc)


install(TARGETS
     ClusterFinder
     ClusterAna_module
     ClusterCheater_module
     Cluster3D_module
     ClusterPCA_module
     DBclusterAna_module
     DumpClusters_module
     EndPointModule_module
     HoughLineFinderAna_module
     HoughLineFinder_module
     LineMerger_module
     SimpleClusterMerger_module
     SmallClusterFinder_module
     fuzzyCluster_module
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


add_subdirectory(RStarTree)
