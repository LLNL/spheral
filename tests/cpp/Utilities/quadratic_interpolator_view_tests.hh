#ifndef __SPHERAL_QUADRATIC_INTERPOLATOR_VIEW_TESTS_HH__
#define __SPHERAL_QUADRATIC_INTERPOLATOR_VIEW_TESTS_HH__

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, DefaultConstructor)
{
  using WORK_EXEC_POLICY = TypeParam;
  Spheral::QuadraticInterpolator q_int;
  auto q_int_v = &q_int;
 
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_EQ(q_int_v->coeffs().size(), 0);
  EXEC_IN_SPACE_END();
}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Initialize)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QuadraticInterpolator q_int;
  q_int.initialize(0,4,{0,1,2});

  auto q_int_v = &q_int;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_EQ(q_int_v->xmin(),          0);
    SPHERAL_ASSERT_EQ(q_int_v->xmax(),          4);
    SPHERAL_ASSERT_EQ(q_int_v->coeffs().size(), 3);
    SPHERAL_ASSERT_EQ(q_int_v->coeffs()[0],     0);
    SPHERAL_ASSERT_EQ(q_int_v->coeffs()[1],     0.5);
    SPHERAL_ASSERT_EQ(q_int_v->coeffs()[2],     0);
  EXEC_IN_SPACE_END();
}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, CopySemantics)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QuadraticInterpolator::CoeffsType c_copy;
  {
    Spheral::QuadraticInterpolator q_int;
    auto q_int_v = &q_int;
    q_int.initialize(0,4,{0,1,2});

    Spheral::QuadraticInterpolator q_int2 = q_int;

    auto q_int2_v = &q_int2;
    auto q_int_v2 = q_int_v;

    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      SPHERAL_ASSERT_EQ(q_int_v->xmin(),          0);
      SPHERAL_ASSERT_EQ(q_int_v->xmax(),          4);
      SPHERAL_ASSERT_EQ(q_int_v->coeffs().size(), 3);
      SPHERAL_ASSERT_EQ(q_int_v->coeffs()[0],     0);
      SPHERAL_ASSERT_EQ(q_int_v->coeffs()[1],     0.5);
      SPHERAL_ASSERT_EQ(q_int_v->coeffs()[2],     0);

      SPHERAL_ASSERT_EQ(q_int_v->xmin(),          q_int2_v->xmin());
      SPHERAL_ASSERT_EQ(q_int_v->xmax(),          q_int2_v->xmax());
      SPHERAL_ASSERT_EQ(q_int_v->coeffs().size(), q_int2_v->coeffs().size());
      SPHERAL_ASSERT_EQ(q_int_v->coeffs()[0],     q_int2_v->coeffs()[0]);
      SPHERAL_ASSERT_EQ(q_int_v->coeffs()[1],     q_int2_v->coeffs()[1]);
      SPHERAL_ASSERT_EQ(q_int_v->coeffs()[2],     q_int2_v->coeffs()[2]);

      SPHERAL_ASSERT_NE(&(q_int_v->coeffs()[0]), &(q_int2_v->coeffs()[0]));
      SPHERAL_ASSERT_EQ(&(q_int_v->coeffs()[0]), &(q_int_v2->coeffs()[0]));
    EXEC_IN_SPACE_END();

    c_copy = deepCopy(q_int.coeffs());
  }
  SPHERAL_ASSERT_EQ(c_copy[0],     0);
  SPHERAL_ASSERT_EQ(c_copy[1],     0.5);
  SPHERAL_ASSERT_EQ(c_copy[2],     0);

  c_copy.free();

}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Equivalence)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QuadraticInterpolator q_int;
  Spheral::QuadraticInterpolator q_int2;
  Spheral::QuadraticInterpolator q_int3;

  auto q_int_v  = &q_int;
  auto q_int2_v = &q_int2;
  auto q_int3_v = &q_int3;

  q_int.initialize(0,4,{0,1,2});
  q_int2.initialize(0,4,{0,1,2});
  q_int3.initialize(0,4,{0,1,3});
 
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_FALSE(q_int_v == q_int2_v);
    SPHERAL_ASSERT_FALSE(q_int_v == q_int3_v);
  EXEC_IN_SPACE_END();
}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, OperatorParen)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QuadraticInterpolator q_int;
  auto q_int_v  = &q_int;

  q_int.initialize(0,4,{0,1,2});
 
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    auto& q = *q_int_v;
    SPHERAL_ASSERT_EQ(q(1), 0.5);
    SPHERAL_ASSERT_EQ(q(2), 1);
    SPHERAL_ASSERT_EQ(q(3), 1.5);
    SPHERAL_ASSERT_EQ(q(4), 2);
   
    SPHERAL_ASSERT_EQ(q(1, 0), 0.5);
    SPHERAL_ASSERT_EQ(q(2, 0), 1);
    SPHERAL_ASSERT_EQ(q(3, 0), 1.5);
    SPHERAL_ASSERT_EQ(q(4, 0), 2);
  EXEC_IN_SPACE_END();
}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Prime)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QuadraticInterpolator q_int;
  auto q_int_v  = &q_int;

  q_int.initialize(0,4,{0,1,2});

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_EQ(q_int_v->prime(1), 0.5);
    SPHERAL_ASSERT_EQ(q_int_v->prime(2), 0.5);
    SPHERAL_ASSERT_EQ(q_int_v->prime(3), 0.5);
    SPHERAL_ASSERT_EQ(q_int_v->prime(4), 0.5);
   
    SPHERAL_ASSERT_EQ(q_int_v->prime(1, 0), 0.5);
    SPHERAL_ASSERT_EQ(q_int_v->prime(2, 0), 0.5);
    SPHERAL_ASSERT_EQ(q_int_v->prime(3, 0), 0.5);
    SPHERAL_ASSERT_EQ(q_int_v->prime(4, 0), 0.5);
  EXEC_IN_SPACE_END();
}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Prime2)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QuadraticInterpolator q_int;
  auto q_int_v  = &q_int;

  q_int.initialize(0,4,{0,1,2});

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_EQ(q_int_v->prime2(1), 0);
    SPHERAL_ASSERT_EQ(q_int_v->prime2(2), 0);
    SPHERAL_ASSERT_EQ(q_int_v->prime2(3), 0);
    SPHERAL_ASSERT_EQ(q_int_v->prime2(4), 0);
   
    SPHERAL_ASSERT_EQ(q_int_v->prime2(1, 0), 0);
    SPHERAL_ASSERT_EQ(q_int_v->prime2(2, 0), 0);
    SPHERAL_ASSERT_EQ(q_int_v->prime2(3, 0), 0);
    SPHERAL_ASSERT_EQ(q_int_v->prime2(4, 0), 0);
  EXEC_IN_SPACE_END();
}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, LowerBound)
{
  using WORK_EXEC_POLICY = TypeParam;

  Spheral::QuadraticInterpolator q_int;
  auto q_int_v  = &q_int;

  q_int.initialize(0,4,{0,1,2});

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_EQ(q_int_v->lowerBound(1), 0);
    SPHERAL_ASSERT_EQ(q_int_v->lowerBound(2), 0);
    SPHERAL_ASSERT_EQ(q_int_v->lowerBound(3), 0);
    SPHERAL_ASSERT_EQ(q_int_v->lowerBound(4), 0);
  EXEC_IN_SPACE_END();
}

#endif //  __SPHERAL_QUADRATIC_INTERPOLATOR_VIEW_TESTS_HH__
