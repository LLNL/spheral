// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

#include "Field/FieldList.hh"
#include "Field/FieldView.hh"
#include "Field/FieldListView.hh"
#include <Utilities/Logger.hh>

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
  GPUCounters fl_count;
  GPUCounters f_count;

  NodeList_t createNodeList(size_t count) {
    return NodeList_t("DataNodeList", count, 0);
  }

  // Increment variables for each action and space
  auto fl_callback() {
    return [&](const chai::PointerRecord *, chai::Action action,
                            chai::ExecutionSpace space) {
    if (action == chai::ACTION_MOVE)
      (space == chai::CPU) ? fl_count.DToHCopies++ : fl_count.HToDCopies++;
    if (action == chai::ACTION_ALLOC)
      (space == chai::CPU) ? fl_count.HNumAlloc++ : fl_count.DNumAlloc++;
    if (action == chai::ACTION_FREE)
      (space == chai::CPU) ? fl_count.HNumFree++ : fl_count.DNumFree++;
    };
  }

  auto f_callback() {
    return [&](const chai::PointerRecord *, chai::Action action,
                            chai::ExecutionSpace space) {
    if (action == chai::ACTION_MOVE)
      (space == chai::CPU) ? f_count.DToHCopies++ : f_count.HToDCopies++;
    if (action == chai::ACTION_ALLOC)
      (space == chai::CPU) ? f_count.HNumAlloc++ : f_count.DNumAlloc++;
    if (action == chai::ACTION_FREE)
      (space == chai::CPU) ? f_count.HNumFree++ : f_count.DNumFree++;
    };
  }
};

// Setting up G Test for FieldList
TYPED_TEST_SUITE_P(FieldListViewTypedTest);
template <typename T> class FieldListViewTypedTest : public FieldListViewTest {};

GPU_TYPED_TEST_P(FieldListViewTypedTest, BasicCapture) {

  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  {
    FieldDouble field("TestField1", nl, val);
    FieldDouble field2("TestField2", nl, val);
    FieldListDouble field_list;

    field_list.appendField(field);
    field_list.appendField(field2);

    const size_t numFields = field_list.size() ;
    auto fl_v = field_list.toView(gpu_this->fl_callback(), gpu_this->f_callback());

    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].size(), N);
      });
  }

  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = {1, 0, 1, 1, 1, 1};
    f_ref_count  = {2, 0, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 1, 0, 1, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

// TODO: Add test for having multiple FL contain the same field

GPU_TYPED_TEST_P(FieldListViewTypedTest, MultiScopeAndTouch) {

  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  {
    FieldDouble field("TestField1", nl, val);
    FieldDouble field2("TestField2", nl, val);
    FieldListDouble field_list;

    field_list.appendField(field);
    field_list.appendField(field2);

    const size_t numFields = field_list.size() ;

    { // Scope 1
    auto fl_v = field_list.toView(gpu_this->fl_callback(), gpu_this->f_callback());

    DEBUG_LOG << "Start Kernel 1";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].size(), N);
      });
    DEBUG_LOG << "Stop Kernel 1";

    fl_v.touch(chai::CPU, true);
    } // Scope 1

    { // Scope 2
    auto fl_v = field_list.toView(gpu_this->fl_callback(), gpu_this->f_callback());

    DEBUG_LOG << "Start Kernel 2";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].size(), N);
      });
    DEBUG_LOG << "Stop Kernel 2";
    } // Scope 2
  }

  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = { 2, 0, 2, 2, 2, 2 };
    f_ref_count  = { 4, 0, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 2, 0, 2, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

GPU_TYPED_TEST_P(FieldListViewTypedTest, MultiScopeNoTouch) {

  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  {
    FieldDouble field("TestField1", nl, val);
    FieldDouble field2("TestField2", nl, val);
    FieldListDouble field_list;

    field_list.appendField(field);
    field_list.appendField(field2);

    const size_t numFields = field_list.size() ;

    { // Scope 1
    auto fl_v = field_list.toView(gpu_this->fl_callback(), gpu_this->f_callback());

    DEBUG_LOG << "Start Kernel 1";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].size(), N);
      });
    DEBUG_LOG << "Stop Kernel 1";
    } // Scope 1

    { // Scope 2
    auto fl_v = field_list.toView(gpu_this->fl_callback(), gpu_this->f_callback());

    DEBUG_LOG << "Start Kernel 2";
    RAJA::forall<TypeParam>(TRS_UINT(0, numFields),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(fl_v.size(), numFields);
        SPHERAL_ASSERT_EQ(fl_v[i].size(), N);
      });
    DEBUG_LOG << "Stop Kernel 2";
    } // Scope 2
  }

  // Counter : { H->D Copy, D->H Copy, H Alloc, D Alloc, H Free, D Free }
  GPUCounters fl_ref_count, f_ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    fl_ref_count = { 2, 0, 2, 2, 2, 2 };
    f_ref_count  = { 2, 0, 0, 2, 0, 2};
  } else {
    fl_ref_count = {0, 0, 2, 0, 2, 0};
    f_ref_count  = {0, 0, 0, 0, 0, 0};
  }

  COMP_COUNTERS(gpu_this->fl_count, fl_ref_count);
  COMP_COUNTERS(gpu_this->f_count,  f_ref_count);
}

// TODO: Add test where 
REGISTER_TYPED_TEST_SUITE_P(FieldListViewTypedTest, BasicCapture, MultiScopeAndTouch, MultiScopeNoTouch);

INSTANTIATE_TYPED_TEST_SUITE_P(FieldListView, FieldListViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
