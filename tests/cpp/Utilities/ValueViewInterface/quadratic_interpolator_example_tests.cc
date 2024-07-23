#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "QInt.hh"

// Setting up G Test for QuadraticInterpolator
template<typename T>
class QIntExampleTypedTest : public::testing::Test {};

// All QuadraticInterpolatorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(QIntExampleTypedTest, EXEC_TYPES);


GPU_TYPED_TEST(QIntExampleTypedTest, SmartCopySemantics)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QInt::CoeffsType c;
  {
    Spheral::QInt qq_int;

    auto qq_int_v = qq_int.toView();
    auto qq_int_v2 = qq_int_v;

    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      SPHERAL_ASSERT_EQ(qq_int_v.xmin(),          0);
      SPHERAL_ASSERT_EQ(qq_int_v.coeffs().size(), 0);
    EXEC_IN_SPACE_END();

    qq_int.initialize(4);
    SPHERAL_ASSERT_EQ(qq_int_v.coeffs().size(), 10);

    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      printf("xmin : %lf\n", qq_int_v.xmin());
      SPHERAL_ASSERT_EQ(qq_int_v.xmin(),          4);
      SPHERAL_ASSERT_EQ(qq_int_v.coeffs().size(), 10);
      SPHERAL_ASSERT_EQ(qq_int_v.coeffs()[9], 0.19);

    EXEC_IN_SPACE_END();

    qq_int.editData(2);

    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      SPHERAL_ASSERT_EQ(qq_int_v.xmin(),          2);
      SPHERAL_ASSERT_EQ(qq_int_v.coeffs()[9], 91);
      printf("xmin : %lf\n", qq_int_v.xmin());
      printf("coeffs[19] : %lf\n", qq_int_v.coeffs()[9]);
    EXEC_IN_SPACE_END();

    for (auto elem : qq_int.coeffs()) std::cout << elem << std::endl;
    c = qq_int.coeffs();
  }
  for (auto elem : c) std::cout << elem << std::endl;
}

