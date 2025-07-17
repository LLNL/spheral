#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"
#include "Neighbor/NodePairList.hh"
#include "Neighbor/NodePairListView.hh"
#include <typeinfo>

using NPIT = Spheral::NodePairIdxType;
using NPLV = Spheral::NodePairListView;
using NPL = Spheral::NodePairList;

static int HToDCopies, DToHCopies, HNumFree, DNumFree, HNumAlloc, DNumAlloc;

static void resetCounters() {
  HToDCopies = 0;
  DToHCopies = 0;
  HNumFree = 0;
  DNumFree = 0;
  HNumAlloc = 0;
  DNumAlloc = 0;
}

// Increment variables for each actionand space
static auto callback = [](const chai::PointerRecord *,
                          chai::Action action,
                          chai::ExecutionSpace space) {
      if (action == chai::ACTION_MOVE) {
        (space == chai::CPU) ? DToHCopies++ : HToDCopies++;
      } else if (action == chai::ACTION_ALLOC) {
        (space == chai::CPU) ? HNumAlloc++ : DNumAlloc++;
      } else if (action == chai::ACTION_FREE) {
        (space == chai::CPU) ? HNumFree++ : DNumFree++;
      }
};

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
  resetCounters();
  const size_t N = 5;
  NPL npl = gpu_this->createContainer(N);
  NPLV npl_v = npl.toView(callback);
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
         npl_v[i].f_couple *= 2.;
       });

  npl_v.move(chai::CPU);

  for (size_t i = 0; i < N; ++i) {
    SPHERAL_ASSERT_EQ(npl[i].i_node, 2*i);
    SPHERAL_ASSERT_EQ(npl[i].i_list, i);
    SPHERAL_ASSERT_EQ(npl[i].f_couple, 2.*(double)i);
  }
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, Touch) {
  resetCounters();
  {
    const size_t N = 5;
    NPL npl = gpu_this->createContainer(N);
    NPLV npl_v = npl.toView(callback);

    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(npl_v.size(), N);
    EXEC_IN_SPACE_END()

    npl[0].i_list = 4; // Modify the data on the host
    npl_v.touch(chai::CPU); // Change the execution space for the MA
    npl_v = npl.toView(callback); // Create a new view

    RAJA::forall<TypeParam>(
         TRS_UINT(0, N),
         [=] SPHERAL_HOST_DEVICE(size_t i) {
           if (i == 0) {
             SPHERAL_ASSERT_EQ(npl_v[i].i_list, 4);
           }
         });
  }
  SPHERAL_ASSERT_EQ(DToHCopies, 0);
  SPHERAL_ASSERT_EQ(HNumAlloc, 0);
  SPHERAL_ASSERT_EQ(HNumFree, 0);
  if (typeid(RAJA::seq_exec) == typeid(TypeParam)) {
    SPHERAL_ASSERT_EQ(HToDCopies, 0);
  } else {
    SPHERAL_ASSERT_EQ(HToDCopies, 2);
    SPHERAL_ASSERT_EQ(DNumAlloc, 1);
    SPHERAL_ASSERT_EQ(DNumFree, 1);
  }
}

// Test constructor from ContainerType, movement to and from device
// and modification on the device
GPU_TYPED_TEST_P(NPLViewTypedTest, Resize) {
  resetCounters();
  {
    const size_t N = 4;
    NPL npl = gpu_this->createContainer(N);
    NPLV npl_v = npl.toView(callback);

    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(npl_v.size(), N);
    EXEC_IN_SPACE_END()

    NPIT nit(4, 4, 4, 4, 4.);
    npl.push_back(nit);
    npl_v = npl.toView(callback);

    RAJA::forall<TypeParam>(
         TRS_UINT(0, npl.size()),
         [=] SPHERAL_HOST_DEVICE(size_t i) {
           SPHERAL_ASSERT_EQ(npl_v.size(), N+1);
           if (i == N) {
             SPHERAL_ASSERT_EQ(npl_v[i].i_list, 4);
             npl_v[i].i_node = 6;
           }
         });
    npl_v.move(chai::CPU);
    SPHERAL_ASSERT_EQ(npl_v[N+1].i_node, npl[N+1].i_node);
  }
  SPHERAL_ASSERT_EQ(HNumAlloc, 0);
  SPHERAL_ASSERT_EQ(HNumFree, 0);
  if (typeid(RAJA::seq_exec) == typeid(TypeParam)) {
    SPHERAL_ASSERT_EQ(HToDCopies, 0);
    SPHERAL_ASSERT_EQ(DToHCopies, 0);
  } else {
    SPHERAL_ASSERT_EQ(DToHCopies, 1);
    SPHERAL_ASSERT_EQ(HToDCopies, 2);
    SPHERAL_ASSERT_EQ(DNumAlloc, 2);
    SPHERAL_ASSERT_EQ(DNumFree, 2);
  }
}

REGISTER_TYPED_TEST_SUITE_P(NPLViewTypedTest, DefaultConstructor,
                            ConstructorFromContainer, Touch, Resize);

INSTANTIATE_TYPED_TEST_SUITE_P(NodePairListView, NPLViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
