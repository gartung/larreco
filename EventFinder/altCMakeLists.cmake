
art_add_module(AggregateEvent_module AggregateEvent_module.cc)

art_add_module(EventCheater_module EventCheater_module.cc)

art_add_module(EventMaker_module EventMaker_module.cc)

install(TARGETS
     AggregateEvent_module
     EventCheater_module
     EventMaker_module
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

