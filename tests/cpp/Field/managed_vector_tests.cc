#include "test-utilities.hh"
#include "test-basic-exec-policies.hh"

#include "Field/SphArray.hh"
#include "chai/managed_ptr.hpp"

using MVDouble = Spheral::ManagedVector<double>;

#define assert_empty_map(IGNORED) ASSERT_EQ(chai::ArrayManager::getInstance()->getPointerMap().size(),0)


// Setting up G Test for ManagedVector
template<typename T>
class ManagedVectorTypedTest : public::testing::Test {};

// All ManagedVectorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(ManagedVectorTypedTest, EXEC_TYPES);


GPU_TYPED_TEST(ManagedVectorTypedTest, DefaultConstructor)
{
  using WORK_EXEC_POLICY = TypeParam;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    MVDouble array;

    SPHERAL_ASSERT_EQ(array.size(), 0);
    SPHERAL_ASSERT_EQ(array.capacity(), 0);
  EXEC_IN_SPACE_END();
}


TEST(ManagedVectorTest, SizeConstructor)
{
  // Size Constructor will allocate initial capacity if elements
  // are below initial_capacity
  MVDouble array(6);
  SPHERAL_ASSERT_EQ(array.size(),     6u);
  SPHERAL_ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);

  // Size Constructor will allocate n elements if elements
  // are above initial_capacity
  MVDouble array2(MVDouble::initial_capacity + 1);
  SPHERAL_ASSERT_EQ(array2.size(),     MVDouble::initial_capacity + 1);
  SPHERAL_ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity * 2);
}



TEST(ManagedVectorTest, IdentityConstructor)
{
  MVDouble array(6, 5);
  SPHERAL_ASSERT_EQ(array.size(),     6u);
  SPHERAL_ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);

  for(auto& elem: array){
    SPHERAL_ASSERT_EQ(elem, 5);
  }
}


GPU_TYPED_TEST(ManagedVectorTypedTest, IdentityConstructor)
{

  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array(6, 5);
  MVDouble array2(6);
  SPHERAL_ASSERT_EQ(array.size(),     6u);
  SPHERAL_ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);

  RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,6),
    [=] RAJA_HOST_DEVICE (unsigned idx){
      array2[idx] = array[idx];
    }
  );
  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,6),
    [=] RAJA_HOST (unsigned idx){
      SPHERAL_ASSERT_EQ(array2[idx], 5);
    }
  );

}


GPU_TYPED_TEST(ManagedVectorTypedTest, CopyConstructor)
{
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array(4);
  MVDouble copy_array(array);
  //MVDouble copy_array = array.slice(0, array.size());

  array.resize(6);

  SPHERAL_ASSERT_EQ(&array[0], &copy_array[0]);
  SPHERAL_ASSERT_EQ(array.capacity(), copy_array.capacity());

  RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,6),
    [=] RAJA_HOST_DEVICE (unsigned i){
      array[i] = i;
    }
  );

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,6),
    [=] RAJA_HOST (unsigned i){
      SPHERAL_ASSERT_EQ(copy_array[i], i);
    }
  );

  SPHERAL_ASSERT_EQ(&array[0], &copy_array[0]);

  array.resize(20);

  SPHERAL_ASSERT_EQ(&array[0], &copy_array[0]);
}

GPU_TYPED_TEST(ManagedVectorTypedTest, AssignmentOperator)
{
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array(6, 5);

  MVDouble copy_array = array;

  RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,6),
    [=] RAJA_HOST_DEVICE (unsigned i){
      array[i] = i;
    }
  );

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,6),
    [=] RAJA_HOST (unsigned i){
      SPHERAL_ASSERT_EQ(copy_array[i], i);
    }
  );

  SPHERAL_ASSERT_EQ(&array[0], &copy_array[0]);
}

GPU_TYPED_TEST(ManagedVectorTypedTest, Equivalence)
{
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array(6, 5);
  MVDouble array2(6, 3);

  MVDouble copy_array = array;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)

    for (size_t i = 0; i < 6u; i++) {
      array[i] = i;
      array2[i] = i;
    }
    SPHERAL_ASSERT_TRUE (copy_array == array);
    SPHERAL_ASSERT_FALSE(copy_array == array2);
  EXEC_IN_SPACE_END()

  EXEC_IN_SPACE_BEGIN(LOOP_EXEC_POLICY)
    ASSERT_TRUE (copy_array == array);
    ASSERT_FALSE(copy_array == array2);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST(ManagedVectorTypedTest, PushBackDefault)
{
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array;
  double val = 5;

  array.push_back(val);

  SPHERAL_ASSERT_EQ(array.size(),     1u);
  SPHERAL_ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
  SPHERAL_ASSERT_EQ(array[0], val);

  MVDouble array2;

  for(size_t i = 0; i < 10; i++) 
    array2.push_back(i);

  for(size_t i = 0; i < array2.size(); i++)
    SPHERAL_ASSERT_EQ(array2[i], i);

  RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0,10),
    [=] RAJA_HOST_DEVICE (unsigned idx){
      SPHERAL_ASSERT_EQ(array2.size(),     10u);
      SPHERAL_ASSERT_EQ(array2[idx], idx);
      array2[idx] *= 2;
    }
  );

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0,6),
    [=] RAJA_HOST (unsigned idx){
      SPHERAL_ASSERT_EQ(array2[idx], idx*2);
    }
  );

  SPHERAL_ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity * 2);
}

TEST(ManagedVectorTest, PushBackMove)
{
  MVDouble array;

  array.push_back(std::move(5));

  SPHERAL_ASSERT_EQ(array.size(),     1u);
  SPHERAL_ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
  SPHERAL_ASSERT_EQ(array[0], 5);

  MVDouble array2;

  for(size_t i = 0; i < 10; i++) 
    array2.push_back(std::move(i));

  for(size_t i = 0; i < array2.size(); i++)
    SPHERAL_ASSERT_EQ(array2[i], i);

  SPHERAL_ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity * 2);
}

GPU_TYPED_TEST(ManagedVectorTypedTest, ResizeLargerNoRealloc)
{
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(array.size(),     0);
  SPHERAL_ASSERT_EQ(array.capacity(), 0);
  EXEC_IN_SPACE_END()

  array.move(chai::CPU);
  array.resize(10);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_EQ(array.size(),     10);
    SPHERAL_ASSERT_EQ(array.capacity(), 16);
  EXEC_IN_SPACE_END()

  MVDouble array2(4);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(array2.size(),     4);
  //SPHERAL_ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity);
  EXEC_IN_SPACE_END()

  array.move(chai::CPU);
  array2.resize(6);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SPHERAL_ASSERT_EQ(array2.size(),     6);
    //SPHERAL_ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST(ManagedVectorTypedTest, ResizeLargerRealloc)
{
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array(4);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(array.size(),     4);
  //SPHERAL_ASSERT_EQ(array.capacity(), 0);
  EXEC_IN_SPACE_END()

  array.move(chai::CPU);
  array.resize(12);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(array.size(),     12);
  SPHERAL_ASSERT_EQ(array.capacity(), 16);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST(ManagedVectorTypedTest, ResizeSmaller)
{
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble array(4);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(array.size(),     4);
  //SPHERAL_ASSERT_EQ(array.capacity(), 0);
  EXEC_IN_SPACE_END()

  array.resize(2);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(array.size(),     2);
  //SPHERAL_ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
  EXEC_IN_SPACE_END()
}

TEST(ManagedVectorTest, Erase)
{
  MVDouble array;
  std::vector<double> check = {0,1,2,3,4,5};

  for (size_t i = 0; i < 6; i++) {
    array.push_back(check[i]);
  }
  SPHERAL_ASSERT_EQ(array.size(), check.size());

  // Erase the last element
  array.erase(array.end() - 1);

  std::vector<double> check2 = {0,1,2,3,4};
  SPHERAL_ASSERT_EQ(array.size(), check2.size());
  for (size_t i = 0; i < array.size(); i++) SPHERAL_ASSERT_EQ(array[i], check2[i]);

  // Erase the 3rd element
  array.erase(array.begin() + 2);

  std::vector<double> check3 = {0,1,3,4};
  SPHERAL_ASSERT_EQ(array.size(), check3.size());
  for (size_t i = 0; i < array.size(); i++) SPHERAL_ASSERT_EQ(array[i], check3[i]);
  
  // Erase the first element
  array.erase(array.begin());

  std::vector<double> check4 = {1,3,4};
  SPHERAL_ASSERT_EQ(array.size(), check4.size());
  for (size_t i = 0; i < array.size(); i++) SPHERAL_ASSERT_EQ(array[i], check4[i]);
}

TEST(ManagedVectorTest, Insert)
{
  MVDouble array;
  std::vector<double> check = {0,1,2,3,4,5};

  for (size_t i = 0; i < 6; i++) {
    array.insert(array.begin() + i, check[i]);
  }
  SPHERAL_ASSERT_EQ(array.size(), check.size());
  for (size_t i = 0; i < array.size(); i++) SPHERAL_ASSERT_EQ(array[i], check[i]);

  // Insert element at end
  array.insert(array.end(), 6);

  std::vector<double> check2 = {0,1,2,3,4,5,6};
  SPHERAL_ASSERT_EQ(array.size(), check2.size());
  for (size_t i = 0; i < array.size(); i++) SPHERAL_ASSERT_EQ(array[i], check2[i]);

  // Erase the 3rd element
  array.insert(array.begin() + 2, 7);

  std::vector<double> check3 = {0,1,7,2,3,4,5,6};
  SPHERAL_ASSERT_EQ(array.size(), check3.size());
  for (size_t i = 0; i < array.size(); i++) SPHERAL_ASSERT_EQ(array[i], check3[i]);
  
  // Erase the first element
  array.insert(array.begin(), -1);

  std::vector<double> check4 = {-1,0,1,7,2,3,4,5,6};
  SPHERAL_ASSERT_EQ(array.size(), check4.size());
  for (size_t i = 0; i < array.size(); i++) SPHERAL_ASSERT_EQ(array[i], check4[i]);
}

TEST(ManagedVectorTest, DeepCopy)
{
  double init_val = 5;
  MVDouble array(6, init_val);

  MVDouble copy_array(array);
  MVDouble deep_copy_array(deepCopy(array));

  for (size_t i = 0; i < 6u; i++) {
    array[i] = 4;
  }

  for (size_t i = 0; i < 6u; i++) {
    ASSERT_NE(array[i], deep_copy_array[i]);
    SPHERAL_ASSERT_EQ(deep_copy_array[i], init_val);
  }

  ASSERT_FALSE(array == deep_copy_array);
  ASSERT_NE(&array[0], &deep_copy_array[0]);
}

