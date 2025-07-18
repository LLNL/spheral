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

static GPUCounters gcounts;

// Increment variables for each action and space
static auto callback = [](const chai::PointerRecord *,
                          chai::Action action,
                          chai::ExecutionSpace space) {
      if (action == chai::ACTION_MOVE) {
        (space == chai::CPU) ?
          gcounts.DToHCopies++ : gcounts.HToDCopies++;
      } else if (action == chai::ACTION_ALLOC) {
        (space == chai::CPU) ?
          gcounts.HNumAlloc++ : gcounts.DNumAlloc++;
      } else if (action == chai::ACTION_FREE) {
        (space == chai::CPU) ?
          gcounts.HNumFree++ : gcounts.DNumFree++;
      }
};

class FieldListViewTest : public ::testing::Test {
public:
  NodeList_t createNodeList(size_t count) {
    return NodeList_t("DataNodeList", count, 0);
  }
};

// Setting up G Test for FieldList
TYPED_TEST_SUITE_P(FieldListViewTypedTest);
template <typename T> class FieldListViewTypedTest : public FieldListViewTest {};

GPU_TYPED_TEST_P(FieldListViewTypedTest, DefaultConstructor) {
  gcounts.resetCounters();
  const int N = 10;
  const double val = 4.;
  std::string field_name = "Field::NodeListValCtor";
  NodeList_t nl = gpu_this->createNodeList(N);
  {
    FieldDouble field(field_name, nl, val);
    FieldDouble field2("Field2", nl, val);
    FieldListDouble field_list;

    field_list.appendField(field);
    field_list.appendField(field2);

    const size_t numFields = field_list.size() ;
    auto fl_v = field_list.toView(callback);

    RAJA::forall<TypeParam>
      (TRS_UINT(0, N),
       [=] SPHERAL_HOST_DEVICE(size_t i) {
         SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
         for (size_t l = 0; l < numFields; ++l) {
           SPHERAL_ASSERT_EQ(fl_v[l].size(), N);
         }
       });
    fl_v.touch(chai::CPU, true);
  }
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count.HToDCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.DNumFree = 1;
  }
  gcounts.compareCounters(ref_count);
}

// TODO: Add test for having multiple FL contain the same field
// TODO: Add test where toView is called in multiple scopes for the same FL
// TODO: Add test where 
REGISTER_TYPED_TEST_SUITE_P(FieldListViewTypedTest, DefaultConstructor);

INSTANTIATE_TYPED_TEST_SUITE_P(FieldListView, FieldListViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
