art_add_module(ShowerCheater_module ShowerCheater_module.cc)

art_add_module(ShowerFinder_module ShowerFinder_module.cc)

art_add_module(ShowerReco_module ShowerReco_module.cc)


install(TARGETS
     ShowerCheater_module
     ShowerFinder_module
     ShowerReco_module
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
