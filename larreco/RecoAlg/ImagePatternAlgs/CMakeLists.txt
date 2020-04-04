# not ready for c2 yet
if( TENSORFLOW_DIR )
  add_subdirectory(Tensorflow)
endif ()
add_subdirectory(DataProvider)
add_subdirectory(job)
add_subdirectory(Keras)
