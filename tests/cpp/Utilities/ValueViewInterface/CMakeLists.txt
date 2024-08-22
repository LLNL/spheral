# This test is based off a basic implementation of Utilities/QuadraticInterpolator.
spheral_add_test(
  NAME value_view_ptr_basic_class_tests
  SOURCES value_view_ptr_basic_class_tests.cc
  #DEBUG_LINKER
)

if(SPHERAL_ENABLE_SCIP)

# This test is based off a basic implementation of Utilities/QuadraticInterpolator.
spheral_add_test(
  NAME value_view_ref_basic_class_tests
  SOURCES value_view_ref_basic_class_tests.cc
  #DEBUG_LINKER
)

# This tests demonstrates CRTP as seen in Kernel.
spheral_add_test(
  NAME value_view_crtp_pattern_tests
  SOURCES value_view_crtp_pattern_tests.cc
)

# This test demonstrates use with a templated abstract base class.
spheral_add_test(
  NAME value_view_template_inheritance_tests
  SOURCES value_view_template_inheritance_tests.cc
)

endif()
