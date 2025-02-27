spheral_add_test(
  NAME quadratic_interpolator_tests
  SOURCES quadratic_interpolator_tests.cc
  DEPENDS_ON Spheral_Utilities
  DEBUG_LINKER
)

spheral_add_test(
  NAME managed_vector_tests
  SOURCES managed_vector_tests.cc
  DEBUG_LINKER
)


add_subdirectory(ValueViewInterface)


