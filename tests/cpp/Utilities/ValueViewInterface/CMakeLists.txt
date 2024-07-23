spheral_add_test(
  NAME quadratic_interpolator_example_tests
  SOURCES quadratic_interpolator_example_tests.cc
  #DEBUG_LINKER
)

spheral_add_test(
  NAME shared_ptr_tests
  SOURCES ManagedSharedPtrTests.cc
  #SOURCES shared_ptr_tests.cc
  #DEBUG_LINKER
)
