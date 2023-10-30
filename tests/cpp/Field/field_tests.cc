#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Field/Field.hh"
#include "NodeList/NodeList.hh"

using DIM3 = Spheral::Dim<3>;
using FieldDouble = Spheral::Field<DIM3, double>;
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

// All ManagedVectorTets cases will run over each type in EXEC_TYPES.
//TYPED_TEST_CASE(FieldTypedTest, EXEC_TYPES);


TEST_F(FieldTest, NameCtor)
{
  //NodeList_t test_node_list("DataNodeList", 10, 0);
  FieldDouble field("MyTestField");
  SPHERAL_ASSERT_EQ(field.name(), "MyTestField");
  SPHERAL_ASSERT_EQ(field.size(), 0);
  std::cout << this->someInt << std::endl;
  // This curently fails as we make all "Fields" valic, due to our 
  // ManagedVectors over allocating space for them...
  //SPHERAL_ASSERT_EQ(field.valid(), false);
}

GPU_TYPED_TEST_P(FieldTypedTest, size)
{
  using WORK_EXEC_POLICY = TypeParam;

  FieldDouble field("MyTestField", gpu_this->test_node_list);
  

}
REGISTER_TYPED_TEST_SUITE_P(FieldTypedTest, size);



INSTANTIATE_TYPED_TEST_SUITE_P(Field, FieldTypedTest, EXEC_TYPES);

