#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"
#include "Neighbor/NodePairList.hh"
#include "Neighbor/NodePairListView.hh"

using NPIT = Spheral::NodePairIdxType;
using NPLV = Spheral::NodePairListView;
using NPL = Spheral::NodePairList;

class NodePairListViewTest : public ::testing::Test {
public:
  // Helper to create a ContainerType with values [start, start+count)
  NPL createContainer(size_t count) {
    NPL npl ;
    for (size_t i = 0; i < count; ++i) {
      NPIT nit(i, i+1, 2*i, 2*i+1, (double)i);
      npl.push_back(nit);
    }
    return npl;
  }
};

// Setting up Typed Test Suite for GeomVector
TYPED_TEST_SUITE_P(NPLViewTypedTest);
template <typename T> class NPLViewTypedTest : public NodePairListViewTest {};

// Test default constructor
GPU_TYPED_TEST_P(NPLViewTypedTest, DefaultConstructor) {
  EXEC_IN_SPACE_BEGIN(TypeParam)
    NPLV npl_v;
    SPHERAL_ASSERT_EQ(npl_v.size(), 0);
  EXEC_IN_SPACE_END()
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, ConstructorFromContainer) {
  const size_t N = 5;
  NPL npl = gpu_this->createContainer(N);
  NPLV npl_v = npl.toView();
  SPHERAL_ASSERT_EQ(npl_v.size(), N);
  SPHERAL_ASSERT_EQ(&(npl_v[0]), &(npl[0]));
  RAJA::forall<TypeParam>(
       TRS_UINT(0, N),
       [=] SPHERAL_HOST_DEVICE(size_t i) {
         NPIT nit = npl_v[i];
         SPHERAL_ASSERT_EQ(nit.i_node, i);
         SPHERAL_ASSERT_EQ(nit.i_list, i+1);
         SPHERAL_ASSERT_EQ(nit.j_node, 2*i);
         SPHERAL_ASSERT_EQ(nit.j_list, 2*i+1);
         SPHERAL_ASSERT_EQ(nit.f_couple, (double)i);
         npl_v[i].i_node *= 2;
         npl_v[i].i_list -= 1;
         npl_v[i].f_couple *= 2.;});
  npl_v.move(chai::CPU);
  for (size_t i = 0; i < N; ++i) {
    SPHERAL_ASSERT_EQ(npl[i].i_node, 2*i);
    SPHERAL_ASSERT_EQ(npl[i].i_list, i);
    SPHERAL_ASSERT_EQ(npl[i].f_couple, 2.*(double)i);
  }
}

REGISTER_TYPED_TEST_SUITE_P(NPLViewTypedTest, DefaultConstructor,
                            ConstructorFromContainer);

INSTANTIATE_TYPED_TEST_SUITE_P(NodePairListView, NPLViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
