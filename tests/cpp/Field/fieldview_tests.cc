// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

#include "Field/FieldView.hh"

/**
 * These are unit tests for Spehral::FieldView with a basic double datatype.
 * Spheral::FieldView is a host/device capable. It is tested using typed
 * tests to check for correct execution on both host and device.
 */

using DIM3 = Spheral::Dim<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using FieldViewDouble = Spheral::FieldView<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;

static GPUCounters gcounts;

// Increment variables for each action and space
static auto callback = [](const chai::PointerRecord *, chai::Action action,
                          chai::ExecutionSpace space) {
  if (action == chai::ACTION_MOVE) {
    if (space == chai::CPU) {
      DEBUG_LOG << "Field Moved to the CPU";
      gcounts.DToHCopies++;
    }
    if (space == chai::GPU) {
      DEBUG_LOG << "Field Moved to the GPU";
      gcounts.HToDCopies++;
    }
  } else if (action == chai::ACTION_ALLOC) {
    if (space == chai::CPU) {
      DEBUG_LOG << "Fieldiew Allocated on the CPU";
      gcounts.HNumAlloc++;
    }
    if (space == chai::GPU) {
      DEBUG_LOG << "Field Allocated on the GPU";
      gcounts.DNumAlloc++;
    }
  } else if (action == chai::ACTION_FREE) {
    if (space == chai::CPU) {
      DEBUG_LOG << "Field DeAllocated on the CPU";
      gcounts.HNumFree++;
    }
    if (space == chai::GPU) {
      DEBUG_LOG << "Field DeAllocated on the GPU";
      gcounts.DNumFree++;
    }
  }
};

class FieldViewTest : public ::testing::Test {
public:
  NodeList_t createNodeList(size_t count) {
    return NodeList_t("DataNodeList", count, 0);
  }
};

// Setting up G Test for Fiedl
TYPED_TEST_SUITE_P(FieldViewTypedTest);
template <typename T> class FieldViewTypedTest : public FieldViewTest {};

/**
 * HOST/Devices CTor test for the Field Ctor that takes a name, nodelist and
 * initial value.
 * - Uses the FieldTest Nodelist constructed for the testing suite above.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, NameNodeListValCtor) {
  gcounts.resetCounters();
  const int N = 10;
  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  int numFields = nl.numFields();
  using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::NodeListValCtor";
    FieldDouble field(field_name, nl, val);
    numFields++;

    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), N);

    auto field_v = field.toView(callback);
    SPHERAL_ASSERT_EQ(field_v.size(), N);

    RAJA::forall<WORK_EXEC_POLICY>
      (TRS_UINT(0, field.size()),
       [=] SPHERAL_HOST_DEVICE(int i) {
         SPHERAL_ASSERT_EQ(field_v[i], val);
       });

    SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
  }
  numFields--;
  SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count.HToDCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.DNumFree = 1;
  }
  gcounts.compareCounters(ref_count);
}

/**
 * Test the multi-view semantics for a copy. If multiple views are made from a
 * single Field then only one copy should be performed as both views will reference
 * the same data.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, MultiViewSemantics) {
  gcounts.resetCounters();
  const int N = 10;
  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  int numFields = nl.numFields();
  using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::MultiViewSemantics";
    FieldDouble field(field_name, nl, val);
    numFields++;

    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), N);

    auto field_v0 = field.toView(callback);
    auto field_v1 = field.toView(callback);
    auto field_v2 = field.toView(callback);
    auto field_v3 = field.toView(callback);
    auto field_v4 = field.toView(callback);
    auto field_v5 = field.toView(callback);
    auto field_v6 = field.toView(callback);
    auto field_v7 = field.toView(callback);
    auto field_v8 = field.toView(callback);
    auto field_v9 = field.toView(callback);
    SPHERAL_ASSERT_EQ(field.size(), N);

    RAJA::forall<WORK_EXEC_POLICY>
      (TRS_UINT(0, field.size()),
       [=] SPHERAL_HOST_DEVICE(int i) {
         SPHERAL_ASSERT_EQ(field_v0[i], val);
         SPHERAL_ASSERT_EQ(field_v1[i], val);
         SPHERAL_ASSERT_EQ(field_v2[i], val);
         SPHERAL_ASSERT_EQ(field_v3[i], val);
         SPHERAL_ASSERT_EQ(field_v4[i], val);
         SPHERAL_ASSERT_EQ(field_v5[i], val);
         SPHERAL_ASSERT_EQ(field_v6[i], val);
         SPHERAL_ASSERT_EQ(field_v7[i], val);
         SPHERAL_ASSERT_EQ(field_v8[i], val);
         SPHERAL_ASSERT_EQ(field_v9[i], val);
       });

    SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
  }
  numFields--;
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count.HToDCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.DNumFree = 1;
  }
  gcounts.compareCounters(ref_count);
  SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
}

/**
 * Resize the field after a copy to the execution space. The Second toView
 * Call should trigger a free of any GPU memory and reassign the FieldView
 * CPU pointer to the underlying vectors new address.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, ResizeField) {
  using WORK_EXEC_POLICY = TypeParam;

  gcounts.resetCounters();
  const int N = 10;
  const double val = 4.;

  NodeList_t nl = gpu_this->createNodeList(N);
  int numFields = nl.numFields();

  {
    std::string field_name = "Field::NodeListValCtor";
    FieldDouble field(field_name, nl, val);
    numFields++;

    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), N);

    auto field_v = field.toView(callback);
    SPHERAL_ASSERT_EQ(field_v.size(), N);

    RAJA::forall<WORK_EXEC_POLICY>
      (TRS_UINT(0, field.size()),
       [=] SPHERAL_HOST_DEVICE(int i) {
         SPHERAL_ASSERT_EQ(field_v[i], val);
       });

    SPHERAL_ASSERT_EQ(field.numInternalElements(), field.numElements());
    nl.numInternalNodes(100);

    field_v = field.toView(callback);
    SPHERAL_ASSERT_EQ(field_v.size(), 100);

    RAJA::forall<WORK_EXEC_POLICY>
      (TRS_UINT(0, field.size()),
       [=] SPHERAL_HOST_DEVICE(int i) {
         if (i < N) {SPHERAL_ASSERT_EQ(field_v[i], val);}
         else { SPHERAL_ASSERT_EQ(field_v[i], 0); }
       });

    SPHERAL_ASSERT_EQ(field.size(), 100);

    SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
  }
  numFields--;

  SPHERAL_ASSERT_EQ(nl.numFields(), numFields);

  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count.HToDCopies = 2;
    ref_count.DNumAlloc = 2;
    ref_count.DNumFree = 2;
  }

  gcounts.compareCounters(ref_count);
}


/**
 * Copy CTor test for the Field.
 * - Test w/ double and GeomPolygon.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, CopyCtor) {
  gcounts.resetCounters();
  const int N = 10;
  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  int numFields = nl.numFields();
  {
    std::string field_name = "Field::CopyCtor";
    FieldDouble field(field_name, nl, val);

    FieldDouble copy_field(field);
    numFields += 2;

    SPHERAL_ASSERT_EQ(copy_field.name(), field_name);
    SPHERAL_ASSERT_EQ(copy_field.size(), N);

    SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);
  }
  numFields -= 2;
  GPUCounters ref_count;
  gcounts.compareCounters(ref_count);
  SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
}

/**
 * Assignment operator of a Field to another Field.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, AssignmentField) {
  gcounts.resetCounters();
  const int N = 10;
  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  int numFields = nl.numFields();
  {
    std::string field_name = "Field::AssignmentField";
    FieldDouble field(field_name, nl, val);

    FieldDouble copy_field("SomeOtherField");
    copy_field = field;
    numFields += 2;

    SPHERAL_ASSERT_EQ(copy_field.size(), N);

    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);

    SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
  }
  numFields -= 2;
  SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
}

/**
 * Assignment operator of a Field to by a std::vector container.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, AssignmentContainerType) {
  gcounts.resetCounters();
  const int N = 10;
  const double val = 4.;
  NodeList_t nl = gpu_this->createNodeList(N);
  int numFields = nl.numFields();
  using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::AssignmentContainer";
    FieldDouble field(field_name, nl, val);
    numFields++;

    using ContainerType = std::vector<double>;
    ContainerType data(N);

    RAJA::forall<LOOP_EXEC_POLICY>
      (TRS_UINT(0, N),
       [&] SPHERAL_HOST(int i) {
         data[i] = i;
       });

    field = data;
    auto field_v = field.toView(callback);

    RAJA::forall<WORK_EXEC_POLICY>
      (TRS_UINT(0, N),
       [=] SPHERAL_HOST_DEVICE(int i) {
         field_v[i] *= 2;
       });

    RAJA::forall<LOOP_EXEC_POLICY>
      (TRS_UINT(0, N),
       [=, &field](int i) {
         SPHERAL_ASSERT_EQ(field_v[i], i * 2);
         SPHERAL_ASSERT_EQ(field[i], i * 2);
       });

    SPHERAL_ASSERT_NE(&field[0], &data[0]);
    SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
  }
  numFields--;
  SPHERAL_ASSERT_EQ(nl.numFields(), numFields);
  GPUCounters ref_count;
  if (typeid(RAJA::seq_exec) != typeid(TypeParam)) {
    ref_count.HToDCopies = 1;
    ref_count.DToHCopies = 1;
    ref_count.DNumAlloc = 1;
    ref_count.DNumFree = 1;
  }
  gcounts.compareCounters(ref_count);
}

REGISTER_TYPED_TEST_SUITE_P(FieldViewTypedTest, NameNodeListValCtor, CopyCtor, MultiViewSemantics,
                            ResizeField, AssignmentField, AssignmentContainerType);

INSTANTIATE_TYPED_TEST_SUITE_P(Field, FieldViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
