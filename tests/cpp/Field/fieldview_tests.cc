#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Field/Field.hh"
#include "Field/FieldView.hh"

using DIM3 = Spheral::Dim<3>;
using FieldDouble     = Spheral::Field<DIM3, double>;
//
using FieldViewDouble = Spheral::FieldView<DIM3, double>;

class FieldViewTest : public::testing::Test{};

// Setting up G Test for Fiedl
TYPED_TEST_SUITE_P(FieldViewTypedTest);
template<typename T>
class FieldViewTypedTest : public FieldViewTest {};

// All ManagedVectorTets cases will run over each type in EXEC_TYPES.
//TYPED_TEST_CASE(FieldViewTypedTest, EXEC_TYPES);


GPU_TYPED_TEST_P(FieldViewTypedTest, DefaultCtor)
{
  using WORK_EXEC_POLICY = TypeParam;

  {

    std::shared_ptr<double> sptr;
    std::cout << sptr.get() << std::endl;

    FieldDouble field("test");
    FieldViewDouble fieldv = field.toView();
    SPHERAL_ASSERT_EQ(fieldv.size(), 0);

    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
      [=] SPHERAL_HOST_DEVICE (int i) {
        SPHERAL_ASSERT_EQ(fieldv.size(), 0);
      });
  }
}


//GPU_TYPED_TEST_P(FieldViewTypedTest, CopyCtor)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  {
//    FieldViewDouble field();
//
//    FieldViewDouble copy_field(field);
//
//    SPHERAL_ASSERT_EQ(copy_field.size(), 10);
//
//    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST_DEVICE (int i) {
//        copy_field[i] == 4;
//      });
//
//    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST_DEVICE (int i) {
//        SPHERAL_ASSERT_EQ(field[i], 4);
//      });
//
//    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 7);
//    SPHERAL_ASSERT_EQ(&field[0], &copy_field[0]);
//  }
//  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
//}


// TODO: ...
//TEST_F(FieldViewTest, AssignmentFieldViewBase)
//{
//  {
//    std::string field_name = "Field::AssignmentFieldBase";
//    FieldDouble field(field_name);
//    SPHERAL_ASSERT_EQ(field.name(), field_name);
//    SPHERAL_ASSERT_EQ(field.size(), 0);
//    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
//  }
//  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
//}

//GPU_TYPED_TEST_P(FieldTypedTest, AssignmentField)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  {
//    std::string field_name = "Field::CopyCtor";
//    FieldDouble field(field_name, gpu_this->test_node_list, 4);
//
//    FieldDouble copy_field = field;
//
//    SPHERAL_ASSERT_EQ(copy_field.name(), "");
//    SPHERAL_ASSERT_EQ(copy_field.size(), 10);
//
//    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST_DEVICE (int i) {
//        copy_field[i] == 4;
//      });
//
//    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST_DEVICE (int i) {
//        SPHERAL_ASSERT_EQ(field[i], 4);
//      });
//
//    SPHERAL_ASSERT_EQ(&field[0], &copy_field[0]);
//
//    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 7);
//  }
//  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
//}
//
//GPU_TYPED_TEST_P(FieldTypedTest, AssignmentContainerType)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  {
//    std::string field_name = "Field::CopyCtor";
//    FieldDouble field(field_name, gpu_this->test_node_list, 4);
//
//    using ContainerType = FieldDouble::ContainerType;
//    ContainerType data(10);
//
//    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST_DEVICE (int i) {
//        data[i] = i;
//      });
//
//    field = data;
//
//    RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST_DEVICE (int i) {
//        field[i] *= 2;
//      });
//
//    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST_DEVICE (int i) {
//        SPHERAL_ASSERT_EQ(field[i], i*2);
//      });
//
//    SPHERAL_ASSERT_EQ(&field[0], &data[0]);
//
//    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
//  }
//  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
//}
//
//TEST_F(FieldTest, AssignmentDataType)
//{
//  {
//    std::string field_name = "Field::NameCtor";
//    FieldDouble field(field_name, this->test_node_list, 4);
//
//    SPHERAL_ASSERT_EQ(field.size(), 10);
//    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
//
//    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST (int i) {
//        SPHERAL_ASSERT_EQ(field[i], 4);
//      });
//
//    field = double(3);
//
//    SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 6);
//    RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,10), 
//      [=] SPHERAL_HOST (int i) {
//        SPHERAL_ASSERT_EQ(field[i], 3);
//      });
//
//
//  }
//  SPHERAL_ASSERT_EQ(this->test_node_list.numFields(), 5);
//}
//
//GPU_TYPED_TEST_P(FieldTypedTest, size)
//{
//  using WORK_EXEC_POLICY = TypeParam;
//
//  {
//    std::string field_name = "Field::size";
//    FieldDouble field(field_name, gpu_this->test_node_list);
//    EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
//      SPHERAL_ASSERT_EQ(field.size(), 10);
//    EXEC_IN_SPACE_END()
//    SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 6);
//  }
//  SPHERAL_ASSERT_EQ(gpu_this->test_node_list.numFields(), 5);
//}

REGISTER_TYPED_TEST_SUITE_P(FieldViewTypedTest,
    DefaultCtor//,
    //CopyCtor,
    ////AssignmentFieldBase,
    //AssignmentField,
    //AssignmentContainerType,
    //size
    );

INSTANTIATE_TYPED_TEST_SUITE_P(FieldView, FieldViewTypedTest, EXEC_TYPES);

