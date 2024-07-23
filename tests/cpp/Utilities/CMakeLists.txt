spheral_add_test(
  NAME quadratic_interpolator_tests
  SOURCES quadratic_interpolator_tests.cc
  DEPENDS_ON Spheral_Utilities
  DEBUG_LINKER
)

add_subdirectory(ValueViewInterface)


