#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

#include "Field/SphArray.hh"

using DIM3 = Spheral::Dim<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using FieldViewDouble = Spheral::FieldView<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;

class FieldTest : public::testing::Test{
  public:
    NodeList_t test_node_list = NodeList_t("DataNodeList", 10, 0);
    int someInt = 5;

    void SetUp() override {

    }

};

// Setting up G Test for Fiedl
TYPED_TEST_SUITE_P(FieldTypedTest);
template<typename T>
class FieldTypedTest : public FieldTest {};

TEST_F(FieldTest, NameCtor)
{
  {
    std::string field_name = "Field::NameCtor";
    FieldDouble field(field_name);
    std::cout << "check\n";
    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 0);
    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

TEST_F(FieldTest, NameNodeListCtor)
{
  {
    std::string field_name = "Field::NodeListCtor";
    FieldDouble field(field_name, this->test_node_list);
    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 10);
    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}



GPU_TYPED_TEST_P(FieldTypedTest, NameNodeListValCtor)
{
  using WORK_EXEC_POLICY = TypeParam;

  {
    std::string field_name = "Field::NodeListValCtor";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);
    SPHERAL_ASSERT_EQ(field.name(), field_name);
    SPHERAL_ASSERT_EQ(field.size(), 10);

    auto field_v = field.toView();

    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        SPHERAL_ASSERT_EQ(field_v[i], 4);
      });
    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}


GPU_TYPED_TEST_P(FieldTypedTest, CopyCtor)
{
  using WORK_EXEC_POLICY = TypeParam;

  {
    std::string field_name = "Field::CopyCtor";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);

    FieldDouble copy_field(field);

    SPHERAL_ASSERT_EQ(copy_field.name(), field_name);
    SPHERAL_ASSERT_EQ(copy_field.size(), 10);

    auto field_v = field.toView();
    auto copy_field_v = copy_field.toView();

    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        SPHERAL_ASSERT_EQ(copy_field_v[i], 4);
        copy_field_v[i] = 5;
      });

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        SPHERAL_ASSERT_EQ(field_v[i], 4);
        SPHERAL_ASSERT_EQ(copy_field_v[i], 5);
      });

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 7);
    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}


GPU_TYPED_TEST_P(FieldTypedTest, FieldView)
{
  using WORK_EXEC_POLICY = TypeParam;

  {
    std::string field_name = "Field::FieldView";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);

    FieldViewDouble field_v(field);

    SPHERAL_ASSERT_EQ(field_v.size(), 10);

    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        field_v[i] = 5;
      });

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        SPHERAL_ASSERT_EQ(field_v[i], 5);
      });

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
    SPHERAL_ASSERT_EQ(&field[0], &field_v[0]);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}


TEST_F(FieldTest, AssignmentFieldBase)
{
  {
    std::string field_name = "Field::AssignmentFieldBase";
    FieldDouble field(field_name, this->test_node_list, 4);
   
    FieldBase* base = &field;

    SPHERAL_ASSERT_EQ(base->name(), field_name);
    //SPHERAL_ASSERT_EQ(field.size(), 0);
    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

GPU_TYPED_TEST_P(FieldTypedTest, AssignmentField)
{
  using WORK_EXEC_POLICY = TypeParam;

  {
    std::string field_name = "Field::AssignmentField";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);

    FieldDouble copy_field = field;

    SPHERAL_ASSERT_EQ(copy_field.name(), field_name);
    SPHERAL_ASSERT_EQ(copy_field.size(), 10);

    auto field_v = field.toView();
    auto copy_field_v = copy_field.toView();

    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        SPHERAL_ASSERT_EQ(copy_field_v[i], 4);
        copy_field_v[i] = 5;
      });

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        SPHERAL_ASSERT_EQ(field_v[i], 4);
        SPHERAL_ASSERT_EQ(copy_field_v[i], 5);
      });

    SPHERAL_ASSERT_NE(&field[0], &copy_field[0]);

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 7);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

GPU_TYPED_TEST_P(FieldTypedTest, AssignmentContainerType)
{
  using WORK_EXEC_POLICY = TypeParam;

  {
    std::string field_name = "Field::AssignmentContainer";
    FieldDouble field(field_name, gpu_this->test_node_list, 4);

    using ContainerType = FieldDouble::ContainerType;
    ContainerType data(10);

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        data[i] = i;
      });

    field = data;
    auto field_v = field.toView();

    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        field_v[i] *= 2;
      });

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        //SPHERAL_ASSERT_EQ(field[i], i); // Should not change we deep copy in lambdas...
        SPHERAL_ASSERT_EQ(field_v[i], i*2);
      });

    SPHERAL_ASSERT_NE(&field[0], &data[0]);

    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

TEST_F(FieldTest, AssignmentDataType)
{
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
  {
    std::string field_name = "Field::AssignmentDataType";
    FieldDouble field(field_name, this->test_node_list, 4);

    SPHERAL_ASSERT_EQ(field.size(), 10);

    auto field_v = field.toView();

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST (int i) {
        SPHERAL_ASSERT_EQ(field_v[i], 4);
      });

    field = double(3);

    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST (int i) {
        SPHERAL_ASSERT_EQ(field_v[i], 3);
      });

  }
  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
}

GPU_TYPED_TEST_P(FieldTypedTest, size)
{
  using WORK_EXEC_POLICY = TypeParam;

  {
    std::string field_name = "Field::size";
    FieldDouble field(field_name, gpu_this->test_node_list);
    auto field_v = field.toView();
    SPHERAL_ASSERT_EQ(field.size(), 10);

    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
      SPHERAL_ASSERT_EQ(field_v.size(), 10);
    EXEC_IN_SPACE_END()
    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
  }
  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
}

REGISTER_TYPED_TEST_SUITE_P(FieldTypedTest,
    NameNodeListValCtor,
    CopyCtor,
    FieldView,
    //AssignmentFieldBase,
    AssignmentField,
    AssignmentContainerType,
    size
    );

INSTANTIATE_TYPED_TEST_SUITE_P(Field, FieldTypedTest, EXEC_TYPES);

