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

class FieldViewTest : public ::testing::Test {
public:
  NodeList_t test_node_list = NodeList_t("DataNodeList", 10, 0);
  int someInt = 5;

  void SetUp() override {}
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
  using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::NodeListValCtor";
    FieldDouble field(field_name, gpu_this->test_node_list, 4.0);

    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 10);

    auto field_v = field.toView();
    SPHERAL_ASSERT_EQ(field.size(), 10);

    RAJA::forall<WORK_EXEC_POLICY>(
        TRS_UINT(0, field.size()),
        [=] SPHERAL_HOST_DEVICE(int i) { SPHERAL_ASSERT_EQ(field_v[i], 4.0); });

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

/**
 * Copy CTor test for the Field.
 * - Test w/ double and GeomPolygon.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, CopyCtor) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::CopyCtor";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);

    FieldDouble copy_field(field);

    SPHERAL_ASSERT_EQ(copy_field.name(), field_name);
    SPHERAL_ASSERT_EQ(copy_field.size(), 10);

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 7);
    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

/**
 * Assignment operator of a Field to another Field.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, AssignmentField) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::AssignmentField";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);

    FieldDouble copy_field("SomeOtherField");
    copy_field = field;

    SPHERAL_ASSERT_EQ(copy_field.size(), 10);

    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 7);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

/**
 * Assignment operator of a Field to by a std::vector container.
 */
GPU_TYPED_TEST_P(FieldViewTypedTest, AssignmentContainerType) {
  using WORK_EXEC_POLICY = TypeParam;
  {
    std::string field_name = "Field::AssignmentContainer";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);

    using ContainerType = std::vector<double>;
    ContainerType data(10);

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 10),
                                   [&] SPHERAL_HOST(int i) { data[i] = i; });

    field = data;
    auto field_v = field.toView();

    RAJA::forall<WORK_EXEC_POLICY>(
        TRS_UINT(0, 10), [=] SPHERAL_HOST_DEVICE(int i) { field_v[i] *= 2; });

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 10), [=, &field](int i) {
      SPHERAL_ASSERT_EQ(field_v[i], i * 2);
      SPHERAL_ASSERT_EQ(field[i], i * 2);
    });

    SPHERAL_ASSERT_NE(&field[0], &data[0]);

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

REGISTER_TYPED_TEST_SUITE_P(FieldViewTypedTest, NameNodeListValCtor, CopyCtor,
                            AssignmentField, AssignmentContainerType);

INSTANTIATE_TYPED_TEST_SUITE_P(Field, FieldViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
