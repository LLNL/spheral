#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Geometry/GeomVector.hh"
#include <Geometry/Box1d.hh>

using Vector = Spheral::GeomVector<3>;

class GeomVectorTest : public ::testing::Test {};

// Setting up Typed Test Suite for GeomVector
TYPED_TEST_SUITE_P(GeomVectorTypedTest);
template <typename T> class GeomVectorTypedTest : public GeomVectorTest {};


GPU_TYPED_TEST_P(GeomVectorTypedTest, Ctor) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(1.0, 2.0, 3.0);
    SPHERAL_ASSERT_EQ(v.x(), 1.0);
    SPHERAL_ASSERT_EQ(v.y(), 2.0);
    SPHERAL_ASSERT_EQ(v.z(), 3.0);
  EXEC_IN_SPACE_END()
}


REGISTER_TYPED_TEST_SUITE_P(GeomVectorTypedTest, Ctor);

INSTANTIATE_TYPED_TEST_SUITE_P(GeomVector, GeomVectorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
