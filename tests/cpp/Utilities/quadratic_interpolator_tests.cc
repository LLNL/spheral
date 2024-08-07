#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Utilities/QuadraticInterpolator.hh"

TEST(QuadraticInterpolatorTest, DefaultConstructor)
{
  Spheral::QuadraticInterpolator q_int;
  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), 0);
}

TEST(QuadraticInterpolatorTest, Initialize)
{
  Spheral::QuadraticInterpolator q_int;
  q_int.initialize(0,4,{0,1,2});

  SPHERAL_ASSERT_EQ(q_int.xmin(),          0);
  SPHERAL_ASSERT_EQ(q_int.xmax(),          4);
  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), 3);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     0);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     0.5);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     0);
}

TEST(QuadraticInterpolatorTest, CopySemantics)
{
  Spheral::QuadraticInterpolator q_int;
  q_int.initialize(0,4,{0,1,2});

  Spheral::QuadraticInterpolator q_int2 = q_int;
  SPHERAL_ASSERT_EQ(q_int.xmin(),          0);
  SPHERAL_ASSERT_EQ(q_int.xmax(),          4);
  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), 3);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     0);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     0.5);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     0);

  SPHERAL_ASSERT_EQ(q_int.xmin(),          q_int2.xmin());
  SPHERAL_ASSERT_EQ(q_int.xmax(),          q_int2.xmax());
  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), q_int2.coeffs().size());
  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     q_int2.coeffs()[0]);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     q_int2.coeffs()[1]);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     q_int2.coeffs()[2]);

  SPHERAL_ASSERT_NE(&(q_int.coeffs()[0]), &(q_int2.coeffs()[0]));
}

TEST(QuadraticInterpolatorTest, AssignmentSemantics)
{
  Spheral::QuadraticInterpolator q_int;
  q_int.initialize(0,4,{0,1,2});

  Spheral::QuadraticInterpolator q_int2;// = q_int;
  q_int2 = q_int;
  SPHERAL_ASSERT_EQ(q_int.xmin(),          0);
  SPHERAL_ASSERT_EQ(q_int.xmax(),          4);
  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), 3);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     0);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     0.5);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     0);

  SPHERAL_ASSERT_EQ(q_int.xmin(),          q_int2.xmin());
  SPHERAL_ASSERT_EQ(q_int.xmax(),          q_int2.xmax());
  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), q_int2.coeffs().size());
  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     q_int2.coeffs()[0]);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     q_int2.coeffs()[1]);
  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     q_int2.coeffs()[2]);

  SPHERAL_ASSERT_NE(&(q_int.coeffs()[0]), &(q_int2.coeffs()[0]));
}

TEST(QuadraticInterpolatorTest, Equivalence)
{
  Spheral::QuadraticInterpolator q_int;
  Spheral::QuadraticInterpolator q_int2;
  Spheral::QuadraticInterpolator q_int3;
  q_int.initialize(0,4,{0,1,2});
  q_int2.initialize(0,4,{0,1,2});
  q_int3.initialize(0,4,{0,1,3});
 
  SPHERAL_ASSERT_TRUE(q_int == q_int2);
  SPHERAL_ASSERT_FALSE(q_int == q_int3);
}

TEST(QuadraticInterpolatorTest, OperatorParen)
{
  Spheral::QuadraticInterpolator q_int;
  q_int.initialize(0,4,{0,1,2});

  SPHERAL_ASSERT_EQ(q_int(1), 0.5);
  SPHERAL_ASSERT_EQ(q_int(2), 1);
  SPHERAL_ASSERT_EQ(q_int(3), 1.5);
  SPHERAL_ASSERT_EQ(q_int(4), 2);
 
  SPHERAL_ASSERT_EQ(q_int(1, 0), 0.5);
  SPHERAL_ASSERT_EQ(q_int(2, 0), 1);
  SPHERAL_ASSERT_EQ(q_int(3, 0), 1.5);
  SPHERAL_ASSERT_EQ(q_int(4, 0), 2);
}
// Setting up G Test for QuadraticInterpolator
template<typename T>
class QuadraticInterpolatorTypedTest : public::testing::Test {};

// All QuadraticInterpolatorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(QuadraticInterpolatorTypedTest, EXEC_TYPES);

//#include "quadratic_interpolator_view_tests.hh"
