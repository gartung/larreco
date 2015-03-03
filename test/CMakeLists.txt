IF(ALT_CMAKE)
INCLUDE(altCMakeLists.cmake)
ELSE()
# various integration tests

include(CetTest)

cet_enable_asserts()

add_subdirectory(HitFinder)
ENDIF()
