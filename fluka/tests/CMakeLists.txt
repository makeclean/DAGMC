message(STATUS "Building FluDAG unit tests")

# Libraries
set(LINK_LIBS)
list(APPEND LINK_LIBS gtest)
list(APPEND LINK_LIBS pthread)
list(APPEND LINK_LIBS fludag)
list(APPEND LINK_LIBS ${DAG_LINK_LIBS})

# Include directories
include_directories(..)
include_directories(${GTEST_INCLUDE_DIR})

# Setup unit tests
set(DRIVER fludag_unit_test_driver.cc)
setup_test(fludag_unit_tests cpp ${DRIVER} "${LINK_LIBS}")
if(NASA)
  setup_test(fludag_nasa_unit_tests cpp ${DRIVER} "${LINK_LIBS}")  
endif(NASA)

# Install h5m files
install(FILES slabs.h5m DESTINATION tests)
if(NASA)
  install(FILES nasa_cancer.h5m DESTINATION tests)
endif(NASA)

# Enable unit tests
enable_testing()
