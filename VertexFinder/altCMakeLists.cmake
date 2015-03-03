art_add_module(AggregateVertexAna_module AggregateVertexAna_module.cc)

art_add_module(AggregateVertex_module AggregateVertex_module.cc)

art_add_module(CornerFinder_module CornerFinder_module.cc)

art_add_module(FeatureVertexFinderAna_module FeatureVertexFinderAna_module.cc)

art_add_module(FeatureVertexFinder_module FeatureVertexFinder_module.cc)

art_add_module(HarrisVertexFinder_module HarrisVertexFinder_module.cc)

art_add_module(PrimaryVertexFinder_module PrimaryVertexFinder_module.cc)

art_add_module(VertexCheater_module VertexCheater_module.cc)

art_add_module(VertexFinder2D_module VertexFinder2D_module.cc)

art_add_module(VertexMatch_module VertexMatch_module.cc)

install(TARGETS
     AggregateVertexAna_module
     AggregateVertex_module
     CornerFinder_module
     FeatureVertexFinderAna_module
     FeatureVertexFinder_module
     HarrisVertexFinder_module
     PrimaryVertexFinder_module
     VertexCheater_module
     VertexFinder2D_module
     VertexMatch_module
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
