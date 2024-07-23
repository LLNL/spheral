spheral_add_test(
  NAME kernel_parallel_inheritance_tests
  SOURCES kernel_parallel_inheritance_tests.cc
)

spheral_add_test(
  NAME field_parallel_inheritance_tests
  SOURCES field_parallel_inheritance_tests.cc
)

spheral_add_test(
  NAME quadratic_interpolator_example_tests
  SOURCES quadratic_interpolator_example_tests.cc
  #DEBUG_LINKER
)
