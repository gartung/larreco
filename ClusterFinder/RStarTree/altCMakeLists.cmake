set(RStarTree_HEADERS
     RStarBoundingBox.h
     RStarTree.h
     RStarVisitor.h
     )

install(FILES ${RStarTree_HEADERS} DESTINATION
     ${CMAKE_INSTALL_INCLUDEDIR}/ClusterFinder/RStarTree 
     COMPONENT Development )
