#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Utilities/QuadraticInterpolator.hh"

// Setting up G Test for QuadraticInterpolator
template<typename T>
class QuadraticInterpolatorTypedTest : public::testing::Test {};

// All QuadraticInterpolatorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(QuadraticInterpolatorTypedTest, EXEC_TYPES);

//TEST(QuadraticInterpolatorTest, DefaultConstructor)
//{
//  Spheral::QuadraticInterpolator q_int;
//  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), 0);
//}
//
//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, DefaultConstructor)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//  Spheral::QuadraticInterpolator q_int;
//  auto q_int_v = q_int.toView();
// 
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs().size(), 0);
//  EXEC_IN_SPACE_END();
//}
//
//
//TEST(QuadraticInterpolatorTest, Initialize)
//{
//  Spheral::QuadraticInterpolator q_int;
//  q_int.initialize(0,4,{0,1,2});
//
//  SPHERAL_ASSERT_EQ(q_int.xmin(),          0);
//  SPHERAL_ASSERT_EQ(q_int.xmax(),          4);
//  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), 3);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     0);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     0.5);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     0);
//}
//
//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Initialize)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  Spheral::QuadraticInterpolator q_int;
//  q_int.initialize(0,4,{0,1,2});
//
//  auto q_int_v = q_int.toView();
//
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_EQ(q_int_v.xmin(),          0);
//    SPHERAL_ASSERT_EQ(q_int_v.xmax(),          4);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs().size(), 3);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs()[0],     0);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs()[1],     0.5);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs()[2],     0);
//  EXEC_IN_SPACE_END();
//}
//
//TEST(QuadraticInterpolatorTest, CopySemantics)
//{
//  Spheral::QuadraticInterpolator q_int;
//  q_int.initialize(0,4,{0,1,2});
//
//  Spheral::QuadraticInterpolator q_int2 = q_int;
//  SPHERAL_ASSERT_EQ(q_int.xmin(),          0);
//  SPHERAL_ASSERT_EQ(q_int.xmax(),          4);
//  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), 3);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     0);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     0.5);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     0);
//
//  SPHERAL_ASSERT_EQ(q_int.xmin(),          q_int2.xmin());
//  SPHERAL_ASSERT_EQ(q_int.xmax(),          q_int2.xmax());
//  SPHERAL_ASSERT_EQ(q_int.coeffs().size(), q_int2.coeffs().size());
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[0],     q_int2.coeffs()[0]);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[1],     q_int2.coeffs()[1]);
//  SPHERAL_ASSERT_EQ(q_int.coeffs()[2],     q_int2.coeffs()[2]);
//
//  SPHERAL_ASSERT_NE(&(q_int.coeffs()[0]), &(q_int2.coeffs()[0]));
//}

GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, SmartCopySemantics)
{
  using WORK_EXEC_POLICY = TypeParam;

  //Spheral::QInt::CoeffsType c;
  {
    Spheral::QInt qq_int;

    auto qq_int_v = qq_int.toView();
    auto qq_int_v2 = qq_int_v;

    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      SPHERAL_ASSERT_EQ(qq_int_v.xmin(),          0);
      SPHERAL_ASSERT_EQ(qq_int_v.coeffs().size(), 0);
    EXEC_IN_SPACE_END();

    qq_int.initialize(4);
    SPHERAL_ASSERT_EQ(qq_int_v.coeffs().size(), 20);

    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      printf("xmin : %lf\n", qq_int_v.xmin());
      //SPHERAL_ASSERT_EQ(qq_int_v.coeffs()[0], 0.1);
      SPHERAL_ASSERT_EQ(qq_int_v.xmin(),          4);
      SPHERAL_ASSERT_EQ(qq_int_v.coeffs().size(), 20);
      SPHERAL_ASSERT_EQ(qq_int_v.coeffs()[19], 0.19);

    EXEC_IN_SPACE_END();


    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      SPHERAL_ASSERT_EQ(qq_int_v.xmin(),          4);
    EXEC_IN_SPACE_END();

    for (auto elem : qq_int.coeffs()) std::cout << elem << std::endl;
    //c = qq_int_v.coeffs();
  }
  //for (auto elem : c) std::cout << elem << std::endl;
}

//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, CopySemantics)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  Spheral::QuadraticInterpolator q_int;
//  auto q_int_v = q_int.toView();
//  q_int.initialize(0,4,{0,1,2});
//
//  Spheral::QInt qq_int;
//  //auto qq_int_v = qq_int.toView();
//
//  //auto q_int_v2 = q_int_v;
//
//  Spheral::QuadraticInterpolator q_int2 = q_int;
//
//  auto q_int2_v = q_int2.toView();
//
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_EQ(q_int_v.xmin(),          0);
//    //SPHERAL_ASSERT_EQ(q_int_v.xmax(),          4);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs().size(), 3);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs()[0],     0);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs()[1],     0.5);
//    SPHERAL_ASSERT_EQ(q_int_v.coeffs()[2],     0);
//
//    //SPHERAL_ASSERT_EQ(q_int_v.xmin(),          q_int2_v.xmin());
//    //SPHERAL_ASSERT_EQ(q_int_v.xmax(),          q_int2_v.xmax());
//    //SPHERAL_ASSERT_EQ(q_int_v.coeffs().size(), q_int2_v.coeffs().size());
//    //SPHERAL_ASSERT_EQ(q_int_v.coeffs()[0],     q_int2_v.coeffs()[0]);
//    //SPHERAL_ASSERT_EQ(q_int_v.coeffs()[1],     q_int2_v.coeffs()[1]);
//    //SPHERAL_ASSERT_EQ(q_int_v.coeffs()[2],     q_int2_v.coeffs()[2]);
//
//    //SPHERAL_ASSERT_NE(&(q_int_v.coeffs()[0]), &(q_int2_v.coeffs()[0]));
//    //SPHERAL_ASSERT_EQ(&(q_int_v.coeffs()[0]), &(q_int_v2.coeffs()[0]));
//  EXEC_IN_SPACE_END();
//}

//TEST(QuadraticInterpolatorTest, Equivalence)
//{
//  Spheral::QuadraticInterpolator q_int;
//  Spheral::QuadraticInterpolator q_int2;
//  Spheral::QuadraticInterpolator q_int3;
//  q_int.initialize(0,4,{0,1,2});
//  q_int2.initialize(0,4,{0,1,2});
//  q_int3.initialize(0,4,{0,1,3});
// 
//  SPHERAL_ASSERT_TRUE(q_int == q_int2);
//  SPHERAL_ASSERT_FALSE(q_int == q_int3);
//}
//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Equivalence)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  Spheral::QuadraticInterpolator q_int;
//  Spheral::QuadraticInterpolator q_int2;
//  Spheral::QuadraticInterpolator q_int3;
//
//  auto q_int_v  = q_int.toView();
//  auto q_int2_v = q_int2.toView();
//  auto q_int3_v = q_int3.toView();
//
//  q_int.initialize(0,4,{0,1,2});
//  q_int2.initialize(0,4,{0,1,2});
//  q_int3.initialize(0,4,{0,1,3});
// 
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_FALSE(q_int_v == q_int2_v);
//    SPHERAL_ASSERT_FALSE(q_int_v == q_int3_v);
//  EXEC_IN_SPACE_END();
//}

//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, OperatorParen)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  Spheral::QuadraticInterpolator q_int;
//  q_int.initialize(0,4,{0,1,2});
//
//  std::cout << q_int(1) << std::endl;
//  std::cout << q_int(2) << std::endl;
//  std::cout << q_int(3) << std::endl;
//  std::cout << q_int(4) << std::endl;
// 
//  std::cout << q_int(1, 0) << std::endl;
//  std::cout << q_int(2, 0) << std::endl;
//  std::cout << q_int(3, 0) << std::endl;
//  std::cout << q_int(4, 0) << std::endl;
// 
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_EQ(q_int(1), 0.5);
//    SPHERAL_ASSERT_EQ(q_int(2), 1);
//    SPHERAL_ASSERT_EQ(q_int(3), 1.5);
//    SPHERAL_ASSERT_EQ(q_int(4), 2);
//   
//    SPHERAL_ASSERT_EQ(q_int(1, 0), 0.5);
//    SPHERAL_ASSERT_EQ(q_int(2, 0), 1);
//    SPHERAL_ASSERT_EQ(q_int(3, 0), 1.5);
//    SPHERAL_ASSERT_EQ(q_int(4, 0), 2);
//  EXEC_IN_SPACE_END();
//}
//
//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Prime)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  Spheral::QuadraticInterpolator q_int;
//  q_int.initialize(0,4,{0,1,2});
//
//  std::cout << q_int.prime(1) << std::endl;
//  std::cout << q_int.prime(2) << std::endl;
//  std::cout << q_int.prime(3) << std::endl;
//  std::cout << q_int.prime(4) << std::endl;
// 
//  std::cout << q_int.prime(1, 0) << std::endl;
//  std::cout << q_int.prime(2, 0) << std::endl;
//  std::cout << q_int.prime(3, 0) << std::endl;
//  std::cout << q_int.prime(4, 0) << std::endl;
// 
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_EQ(q_int.prime(1), 0.5);
//    SPHERAL_ASSERT_EQ(q_int.prime(2), 0.5);
//    SPHERAL_ASSERT_EQ(q_int.prime(3), 0.5);
//    SPHERAL_ASSERT_EQ(q_int.prime(4), 0.5);
//   
//    SPHERAL_ASSERT_EQ(q_int.prime(1, 0), 0.5);
//    SPHERAL_ASSERT_EQ(q_int.prime(2, 0), 0.5);
//    SPHERAL_ASSERT_EQ(q_int.prime(3, 0), 0.5);
//    SPHERAL_ASSERT_EQ(q_int.prime(4, 0), 0.5);
//  EXEC_IN_SPACE_END();
//}
//
//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, Prime2)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  Spheral::QuadraticInterpolator q_int;
//  q_int.initialize(0,4,{0,1,2});
//
//  std::cout << q_int.prime2(1) << std::endl;
//  std::cout << q_int.prime2(2) << std::endl;
//  std::cout << q_int.prime2(3) << std::endl;
//  std::cout << q_int.prime2(4) << std::endl;
// 
//  std::cout << q_int.prime2(1, 0) << std::endl;
//  std::cout << q_int.prime2(2, 0) << std::endl;
//  std::cout << q_int.prime2(3, 0) << std::endl;
//  std::cout << q_int.prime2(4, 0) << std::endl;
// 
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_EQ(q_int.prime2(1), 0);
//    SPHERAL_ASSERT_EQ(q_int.prime2(2), 0);
//    SPHERAL_ASSERT_EQ(q_int.prime2(3), 0);
//    SPHERAL_ASSERT_EQ(q_int.prime2(4), 0);
//   
//    SPHERAL_ASSERT_EQ(q_int.prime2(1, 0), 0);
//    SPHERAL_ASSERT_EQ(q_int.prime2(2, 0), 0);
//    SPHERAL_ASSERT_EQ(q_int.prime2(3, 0), 0);
//    SPHERAL_ASSERT_EQ(q_int.prime2(4, 0), 0);
//  EXEC_IN_SPACE_END();
//}
//
//GPU_TYPED_TEST(QuadraticInterpolatorTypedTest, LowerBound)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  Spheral::QuadraticInterpolator q_int;
//  q_int.initialize(0,4,{0,1,2});
//
//  std::cout << q_int.lowerBound(1) << std::endl;
//  std::cout << q_int.lowerBound(2) << std::endl;
//  std::cout << q_int.lowerBound(3) << std::endl;
//  std::cout << q_int.lowerBound(4) << std::endl;
// 
//  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//    SPHERAL_ASSERT_EQ(q_int.lowerBound(1), 0);
//    SPHERAL_ASSERT_EQ(q_int.lowerBound(2), 0);
//    SPHERAL_ASSERT_EQ(q_int.lowerBound(3), 0);
//    SPHERAL_ASSERT_EQ(q_int.lowerBound(4), 0);
//  EXEC_IN_SPACE_END();
//}
