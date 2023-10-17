#include "gtest/gtest.h"


#include "Field/SphArray.hh"

using MVDouble = Spheral::ManagedVector<double>;


TEST(ManagedVector, DefaultConstructor)
{
  MVDouble array;

  // Size should remain 0
  ASSERT_EQ(array.size(),     0u);
  // This should allocate the reserved capacity value.
  ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
}


TEST(ManagedVector, SizeConstructor)
{
  // Size Constructor will allocate initial capacity if elements
  // are below initial_capacity
  MVDouble array(6);
  ASSERT_EQ(array.size(),     6u);
  ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);

  // Size Constructor will allocate n elements if elements
  // are above initial_capacity
  MVDouble array2(MVDouble::initial_capacity + 1);
  ASSERT_EQ(array2.size(),     MVDouble::initial_capacity + 1);
  ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity + 1);
}



TEST(ManagedVector, IdentityConstructor)
{
  MVDouble array(6, 5);
  ASSERT_EQ(array.size(),     6u);
  ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);

  for(auto& elem: array){
    ASSERT_EQ(elem, 5);
  }
}

TEST(ManagedVector, CopyConstructor)
{
  MVDouble array(6, 5);

  MVDouble copy_array(array);

  for (size_t i = 0; i < 6u; i++) {
    array[i] = i;
  }

  for (size_t i = 0; i < 6u; i++) {
    ASSERT_EQ(copy_array[i], i);
  }

  ASSERT_EQ(&array[0], &copy_array[0]);
}

TEST(ManagedVector, AssignmentOperator)
{
  MVDouble array(6, 5);

  MVDouble copy_array = array;

  for (size_t i = 0; i < 6u; i++) {
    array[i] = i;
  }

  for (size_t i = 0; i < 6u; i++) {
    ASSERT_EQ(copy_array[i], i);
  }

  ASSERT_EQ(&array[0], &copy_array[0]);
}

TEST(ManagedVector, Equivalence)
{
  MVDouble array(6, 5);
  MVDouble array2(6, 5);

  MVDouble copy_array = array;

  for (size_t i = 0; i < 6u; i++) {
    array[i] = i;
    array2[i] = i;
  }

  ASSERT_TRUE (copy_array == array);
  ASSERT_FALSE(copy_array == array2);
}

TEST(ManagedVector, PushBackDefault)
{
  MVDouble array;
  double val = 5;

  array.push_back(val);

  ASSERT_EQ(array.size(),     1u);
  ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
  ASSERT_EQ(array[0], val);

  MVDouble array2;

  for(size_t i = 0; i < 10; i++) 
    array2.push_back(i);

  for(size_t i = 0; i < array2.size(); i++)
    ASSERT_EQ(array2[i], i);

  ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity + (MVDouble::initial_capacity / 2) );
}

TEST(ManagedVector, PushBackMove)
{
  MVDouble array;

  array.push_back(std::move(5));

  ASSERT_EQ(array.size(),     1u);
  ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
  ASSERT_EQ(array[0], 5);

  MVDouble array2;

  for(size_t i = 0; i < 10; i++) 
    array2.push_back(std::move(i));

  for(size_t i = 0; i < array2.size(); i++)
    ASSERT_EQ(array2[i], i);

  ASSERT_EQ(array2.capacity(), MVDouble::initial_capacity + (MVDouble::initial_capacity / 2) );
}

TEST(ManagedVector, ResizeLargerNoRealloc)
{
  MVDouble array(4);

  array.resize(6);

  ASSERT_EQ(array.size(),     6);
  ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
}

TEST(ManagedVector, ResizeLargerRealloc)
{
  MVDouble array(4);

  array.resize(12);

  ASSERT_EQ(array.size(),     12);
  ASSERT_EQ(array.capacity(), 12);
}

TEST(ManagedVector, ResizeSmaller)
{
  MVDouble array(4);

  array.resize(2);

  ASSERT_EQ(array.size(),     2);
  ASSERT_EQ(array.capacity(), MVDouble::initial_capacity);
}

TEST(ManagedVector, Erase)
{
  MVDouble array;
  std::vector<double> check = {0,1,2,3,4,5};

  for (size_t i = 0; i < 6; i++) {
    array.push_back(check[i]);
  }
  ASSERT_EQ(array.size(), check.size());

  // Erase the last element
  array.erase(array.end() - 1);

  std::vector<double> check2 = {0,1,2,3,4};
  ASSERT_EQ(array.size(), check2.size());
  for (size_t i = 0; i < array.size(); i++) ASSERT_EQ(array[i], check2[i]);

  // Erase the 3rd element
  array.erase(array.begin() + 2);

  std::vector<double> check3 = {0,1,3,4};
  ASSERT_EQ(array.size(), check3.size());
  for (size_t i = 0; i < array.size(); i++) ASSERT_EQ(array[i], check3[i]);
  
  // Erase the first element
  array.erase(array.begin());

  std::vector<double> check4 = {1,3,4};
  ASSERT_EQ(array.size(), check4.size());
  for (size_t i = 0; i < array.size(); i++) ASSERT_EQ(array[i], check4[i]);
}

TEST(ManagedVector, DeepCopy)
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
    ASSERT_EQ(deep_copy_array[i], init_val);
  }

  ASSERT_FALSE(array == deep_copy_array);
  ASSERT_NE(&array[0], &deep_copy_array[0]);
}

