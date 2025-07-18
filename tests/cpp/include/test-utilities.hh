#ifndef SPHERAL_TEST_UTILITIES_HH
#define SPHERAL_TEST_UTILITIES_HH

#include "RAJA/RAJA.hpp"
#include "assert.h"
#include "config.hh"
#include "gtest/gtest.h"

using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;
using LOOP_EXEC_POLICY = RAJA::seq_exec;

#define EXEC_IN_SPACE_BEGIN(POL)                                               \
  RAJA::forall<POL>(TRS_UINT(0,1), [=] SPHERAL_HOST_DEVICE (int) {

#define EXEC_IN_SPACE_END()                                                    \
  });

#if defined(SPHERAL_ENABLE_HIP) && defined(__HIPCC__)
#define SPHERAL_GPU_ACTIVE
#endif // SPHERAL_ENABLE_CUDA && __CUDACC__

// Cannot be called on device
#define SPHERAL_ASSERT_EQ_MSG(LHS, RHS) ASSERT_EQ(LHS, RHS);
#if !defined(SPHERAL_GPU_ACTIVE)

#define SPHERAL_ASSERT_EQ(LHS, RHS) ASSERT_EQ(LHS, RHS);
#define SPHERAL_ASSERT_EQ_MSG(LHS, RHS, MSG) ASSERT_EQ(LHS, RHS);
#define SPHERAL_ASSERT_NE(LHS, RHS) ASSERT_NE(LHS, RHS);
#define SPHERAL_ASSERT_TRUE(VAL) ASSERT_TRUE(VAL);
#define SPHERAL_ASSERT_FALSE(VAL) ASSERT_FALSE(VAL);
#define SPHERAL_ASSERT_FLOAT_EQ(LHS, RHS) ASSERT_FLOAT_EQ(LHS, RHS);

#else

#define SPHERAL_ASSERT_EQ_MSG(LHS, RHS, MSG)                                   \
  if (LHS != RHS) {                                                            \
    printf(MSG);                                                               \
    assert(0);                                                                 \
  }

#define SPHERAL_ASSERT_EQ(LHS, RHS)                                            \
  if (LHS != RHS) {                                                            \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#define SPHERAL_ASSERT_NE(LHS, RHS)                                            \
  if (LHS == RHS) {                                                            \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#define SPHERAL_ASSERT_TRUE(VAL)                                               \
  if (!(VAL)) {                                                                \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#define SPHERAL_ASSERT_FALSE(VAL)                                              \
  if ((VAL)) {                                                                 \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#define SPHERAL_ASSERT_FLOAT_EQ(LHS, RHS)                                      \
  if (fabs(LHS - RHS) > 1e-5) {                                                \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#endif

// Macro used throughout LLNLProjects to get around calling
// HOST_DEVICE lambdas from within the "Testing" class directly
#define GPU_TYPED_TEST(X, Y)                                                   \
  template <typename TypeParam> static void gpu_test_##X##Y();                 \
  TYPED_TEST(X, Y) { gpu_test_##X##Y<TypeParam>(); }                           \
  template <typename TypeParam> static void gpu_test_##X##Y()

#define GPU_TYPED_TEST_P(X, Y)                                                 \
  template <typename TypeParam, typename TestFixture>                          \
  static void gpu_test_##X##Y(TestFixture *gpu_this);                          \
  TYPED_TEST_P(X, Y) { gpu_test_##X##Y<TypeParam>(this); }                     \
  template <typename TypeParam, typename TestFixture>                          \
  static void gpu_test_##X##Y(TestFixture *gpu_this)

struct GPUCounters {
  int HToDCopies, DToHCopies;
  int HNumAlloc, DNumAlloc;
  int HNumFree, DNumFree;

  // Constructor to initialize all counters to zero
  GPUCounters() {
    resetCounters();
  }

  void resetCounters() {
    HToDCopies = 0, DToHCopies = 0, HNumAlloc = 0;
    DNumAlloc = 0, HNumFree = 0, DNumFree = 0;
  }

  // Function to compare with another GPUCounter struct
  void compareCounters(const GPUCounters& ref) {
    SPHERAL_ASSERT_EQ_MSG(HToDCopies, ref.HToDCopies);
    SPHERAL_ASSERT_EQ_MSG(DToHCopies, ref.DToHCopies);
    SPHERAL_ASSERT_EQ_MSG(HNumAlloc, ref.HNumAlloc);
    SPHERAL_ASSERT_EQ_MSG(DNumAlloc, ref.DNumAlloc);
    SPHERAL_ASSERT_EQ_MSG(HNumFree, ref.HNumFree);
    SPHERAL_ASSERT_EQ_MSG(DNumFree, ref.DNumFree);
  }
};

#endif // SPHERAL_TEST_UTILITIES_HH
