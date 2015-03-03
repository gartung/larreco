art_add_module(TTSpacePointFinder_module TTSpacePointFinder_module.cc)

install(TARGETS
     TTSpacePointFinder_module
     EXPORT larrecoLibraries
     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
     COMPONENT Runtime 
     )

install(FILES  spptfindermodules.fcl DESTINATION job COMPONENT Runtime)
