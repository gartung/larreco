INSTALL(FILES
  ArtCPack.cmake
  ArtDictionary.cmake
  ArtMake.cmake
  BuildDictionary.cmake
  BuildPlugins.cmake
  CetParseArgs.cmake
  CetRegexEscape.cmake
  CetRootCint.cmake
  CetTest.cmake
  CheckClassVersion.cmake
  FindCppUnit.cmake
  FindSQLite3.cmake
  FindTBB.cmake
  portableArtDictionary.cmake
  portableCheckClassVersion.cmake
  artInternalTools.cmake
  artTools.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/art-${art_VERSION}
  COMPONENT Development
  )
