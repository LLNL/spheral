// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Utilities/CubicHermiteInterpolator.hh"
#include <Utilities/Logger.hh>
#include <functional>

using QI = Spheral::CubicHermiteInterpolator;

class CubicHermiteInterpolatorTest : public ::testing::Test {
public:
  const double xmin = 10.;
  const double xmax = 100.;
  SPHERAL_HOST_DEVICE static double func(const double x) {
    return 1. + x*(2. + 3.*x);
  }
  std::vector<double> makeVec(size_t N) {
    double xstep = (xmax - xmin)/((double)N - 1.);
    std::vector<double> yvals(N);
    for (size_t i = 0; i < N; ++i) {
      double x = xmin + xstep*(double)i;
      yvals[i] = func(x);
    }
    return yvals;
  }
};

// Setting up G Test for FieldList
TYPED_TEST_SUITE_P(CubicHermiteInterpolatorTypedTest);
template <typename T> class CubicHermiteInterpolatorTypedTest : public CubicHermiteInterpolatorTest {};

// Test multiple FieldLists holding the same Field
GPU_TYPED_TEST_P(CubicHermiteInterpolatorTypedTest, FuncCtorTest) {
  const size_t NV = 41;
  const double xmin = gpu_this->xmin;
  const double xmax = gpu_this->xmax;
  QI qih(xmin, xmax, NV, gpu_this->func);
  {
    size_t N = qih.size();
    Spheral::CHIBase qi = qih.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(qi.size(), N);
    EXEC_IN_SPACE_END()
    const double xstep = (xmax - xmin)/((double)NV - 1.);
    RAJA::forall<TypeParam>(TRS_UINT(0, NV),
      [=] (size_t i) {
        double x = xmin + xstep*(double)i;
        double rval = gpu_this->func(x);
        double ival = qi(x);
        SPHERAL_ASSERT_FLOAT_EQ(rval, ival);
      });
  }
}

GPU_TYPED_TEST_P(CubicHermiteInterpolatorTypedTest, VecCtorTest) {
  const size_t NV = 41;
  std::vector<double> yvals = gpu_this->makeVec(NV);
  const double xmin = gpu_this->xmin;
  const double xmax = gpu_this->xmax;
  QI qih(xmin, xmax, yvals);
  size_t N = qih.size();
  Spheral::CHIBase qi = qih.view();
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_EQ(qi.size(), N);
  EXEC_IN_SPACE_END()
  const double xstep = (xmax - xmin)/((double)NV - 1.);
  RAJA::forall<TypeParam>(TRS_UINT(0, NV),
    [=] (size_t i) {
      double x = xmin + xstep*(double)i;
      double rval = gpu_this->func(x);
      double ival = qi(x);
      SPHERAL_ASSERT_FLOAT_EQ(rval, ival);
    });
}

REGISTER_TYPED_TEST_SUITE_P(CubicHermiteInterpolatorTypedTest, FuncCtorTest,
                            VecCtorTest);

INSTANTIATE_TYPED_TEST_SUITE_P(CubicHermiteInterpolator, CubicHermiteInterpolatorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
