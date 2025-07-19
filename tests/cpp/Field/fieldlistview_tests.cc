// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

#include "Field/FieldList.hh"
#include "Field/FieldView.hh"
#include "Field/FieldListView.hh"

using DIM3 = Spheral::Dim<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using FieldViewDouble = Spheral::FieldView<DIM3, double>;
using FieldListDouble = Spheral::FieldList<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;

// Default Testing Size.
static constexpr int N = 100;

class FieldListViewTest : public ::testing::Test {
public:
  GPUCounters gcounts;

  NodeList_t createNodeList(size_t count) {
    return NodeList_t("DataNodeList", count, 0);
  }

  // Increment variables for each action and space
  auto callback() {
    return [&](const chai::PointerRecord *, chai::Action action,
                            chai::ExecutionSpace space) {
    if (action == chai::ACTION_MOVE)
      (space == chai::CPU) ? gcounts.DToHCopies++ : gcounts.HToDCopies++;
    if (action == chai::ACTION_ALLOC)
      (space == chai::CPU) ? gcounts.HNumAlloc++ : gcounts.DNumAlloc++;
    if (action == chai::ACTION_FREE)
      (space == chai::CPU) ? gcounts.HNumFree++ : gcounts.DNumFree++;
    };
  }
};

// Setting up G Test for FieldList
TYPED_TEST_SUITE_P(FieldListViewTypedTest);
template <typename T> class FieldListViewTypedTest : public FieldListViewTest {};

GPU_TYPED_TEST_P(FieldListViewTypedTest, DefaultConstructor) {

  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  {
    FieldDouble field("TestField1", nl, val);
    FieldDouble field2("TestField2", nl, val);
    FieldListDouble field_list;

    field_list.appendField(field);
    field_list.appendField(field2);

    const size_t numFields = field_list.size() ;
    auto fl_v = field_list.toView(gpu_this->callback());

    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].size(), N);
      });
    fl_v.touch(chai::CPU, true);
  }

  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count.HToDCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.HNumFree = 1;
    ref_count.DNumFree = 1;
  } else {
    ref_count.HNumFree = 1;
  }
  COMP_COUNTERS(gpu_this->gcounts, ref_count);
}

// TODO: Add test for having multiple FL contain the same field
// TODO: Add test where toView is called in multiple scopes for the same FL
// TODO: Add test where 
REGISTER_TYPED_TEST_SUITE_P(FieldListViewTypedTest, DefaultConstructor);

INSTANTIATE_TYPED_TEST_SUITE_P(FieldListView, FieldListViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
