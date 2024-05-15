#ifndef SPHERAL_TEST_UTIILITIES_HH
#define SPHERAL_TEST_UTIILITIES_HH


#include "gtest/gtest.h"
#include "assert.h"
#include "RAJA/RAJA.hpp"
#include "config.hh"

using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;
using LOOP_EXEC_POLICY = RAJA::seq_exec;

#define EXEC_IN_SPACE_BEGIN(POL) \
  RAJA::forall<POL>(TRS_UINT(0,1), [=] SPHERAL_HOST_DEVICE (int) {

#define EXEC_IN_SPACE_END() \
  });


template<typename T, typename U>
inline
SPHERAL_HOST_DEVICE
void SPHERAL_ASSERT_EQ(T const& LHS, U const& RHS) {
#if !defined(SPHERAL_GPU_ACTIVE)
  ASSERT_EQ(LHS, RHS);
#else
  if (LHS != RHS) {
    printf("ERROR @ cuda_assert\n");
    assert(0); 
  }
#endif
}

template<typename T, typename U>
SPHERAL_HOST_DEVICE
void SPHERAL_ASSERT_NE(T const& LHS, U const& RHS) {
#if !defined(SPHERAL_GPU_ACTIVE)
  ASSERT_NE(LHS, RHS);
#else
  if (LHS == RHS) {
    printf("ERROR @ cuda_assert\n");
    assert(0); 
  }
#endif
}

SPHERAL_HOST_DEVICE
void SPHERAL_ASSERT_TRUE(bool result) {
#if !defined(SPHERAL_GPU_ACTIVE)
  ASSERT_TRUE(result);
#else
  if (!result) {
    printf("ERROR @ cuda_assert\n");
    assert(0); 
  }
#endif
}

SPHERAL_HOST_DEVICE
void SPHERAL_ASSERT_FALSE(bool result) {
#if !defined(SPHERAL_GPU_ACTIVE)
  ASSERT_FALSE(result);
#else
  if (result) {
    printf("ERROR @ cuda_assert\n");
    assert(0); 
  }
#endif
}


// Macro used throughout LLNLProjects to get around calling
// HOST_DEVICE lambdas from within the "Testing" class directly
#define GPU_TYPED_TEST(X, Y)              \
  template<typename TypeParam> \
  static void gpu_test_##X##Y();    \
  TYPED_TEST(X, Y) { gpu_test_##X##Y<TypeParam>(); } \
  template<typename TypeParam> \
  static void gpu_test_##X##Y()

#define GPU_TYPED_TEST_P(X, Y)              \
  template<typename TypeParam, typename TestFixture> \
  static void gpu_test_##X##Y(TestFixture* gpu_this);    \
  TYPED_TEST_P(X, Y) { gpu_test_##X##Y<TypeParam>(this); } \
  template<typename TypeParam, typename TestFixture> \
  static void gpu_test_##X##Y(TestFixture* gpu_this)


#endif // SPHERAL_TEST_UTIILITIES_HH
