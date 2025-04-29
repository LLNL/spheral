#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Utilities/ManagedVector.hh"
#include "Geometry/GeomPolygon.hh"

using MVDouble = Spheral::ManagedVector<double>;
using MVValType = typename MVDouble::value_type;

#define assert_empty_map(IGNORED)                                              \
  ASSERT_EQ(chai::ArrayManager::getInstance()->getPointerMap().size(), 0)

// Setting up G Test for ManagedVector
template <typename T> class ManagedVectorTypedTest : public ::testing::Test {};

// All ManagedVectorTets cases will run over each type in EXEC_TYPES.
TYPED_TEST_CASE(ManagedVectorTypedTest, EXEC_TYPES);

GPU_TYPED_TEST(ManagedVectorTypedTest, DefaultConstructor) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr;

  SPHERAL_ASSERT_EQ(arr.size(), 0);
  SPHERAL_ASSERT_EQ(arr.capacity(), 0);

  arr.free();
}

TEST(ManagedVectorTest, SizeConstructor) {
  // Size Constructor will allocate initial capacity if elements
  // are below initial_capacity
  MVDouble arr(6);
  SPHERAL_ASSERT_EQ(arr.size(), 6u);
  SPHERAL_ASSERT_EQ(arr.capacity(), MVDouble::initial_capacity);

  // Size Constructor will allocate n elements if elements
  // are above initial_capacity
  MVDouble arr2(MVDouble::initial_capacity + 1);
  SPHERAL_ASSERT_EQ(arr2.size(), MVDouble::initial_capacity + 1);
  SPHERAL_ASSERT_EQ(arr2.capacity(), MVDouble::initial_capacity * 2);

  arr.free();
  arr2.free();
}

TEST(ManagedVectorTest, IdentityConstructor) {
  MVDouble arr(6, 5);
  SPHERAL_ASSERT_EQ(arr.size(), 6u);
  SPHERAL_ASSERT_EQ(arr.capacity(), MVDouble::initial_capacity);

  for (auto &elem : arr) {
    SPHERAL_ASSERT_EQ(elem, 5);
  }
  arr.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, IdentityConstructor) {

  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr(6, 5);
  MVDouble arr2(6);
  SPHERAL_ASSERT_EQ(arr.size(), 6u);
  SPHERAL_ASSERT_EQ(arr.capacity(), MVDouble::initial_capacity);

  RAJA::forall<WORK_EXEC_POLICY>(
      TRS_UINT(0, 6),
      [=] RAJA_HOST_DEVICE(unsigned idx) { arr2[idx] = arr[idx]; });
  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 6), [=] RAJA_HOST(unsigned idx) {
    SPHERAL_ASSERT_EQ(arr2[idx], 5);
  });

  arr.free();
  arr2.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, CopyConstructor) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr(4);
  MVDouble copy_arr(arr);
  // MVDouble copy_arr = arr.slice(0, arr.size());

  arr.resize(6);

  SPHERAL_ASSERT_EQ(&arr[0], &copy_arr[0]);
  SPHERAL_ASSERT_EQ(arr.capacity(), copy_arr.capacity());

  RAJA::forall<WORK_EXEC_POLICY>(
      TRS_UINT(0, 6), [=] RAJA_HOST_DEVICE(unsigned i) { arr[i] = i; });

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 6), [=] RAJA_HOST(unsigned i) {
    SPHERAL_ASSERT_EQ(copy_arr[i], i);
  });

  SPHERAL_ASSERT_EQ(&arr[0], &copy_arr[0]);

  arr.resize(20);

  SPHERAL_ASSERT_EQ(&arr[0], &copy_arr[0]);

  copy_arr.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, DeepCopy) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr(6);

  // MVDouble copy_arr = arr.slice(0, arr.size());

  RAJA::forall<WORK_EXEC_POLICY>(
      TRS_UINT(0, arr.size()),
      [=] RAJA_HOST_DEVICE(unsigned i) { arr[i] = i; });

  MVDouble copy_arr = arr.clone();

  SPHERAL_ASSERT_NE(&arr[0], &copy_arr[0]);

  SPHERAL_ASSERT_EQ(arr.size(), copy_arr.size());
  SPHERAL_ASSERT_EQ(arr.capacity(), copy_arr.capacity());

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 6), [=] RAJA_HOST(unsigned i) {
    SPHERAL_ASSERT_EQ(copy_arr[i], i);
    SPHERAL_ASSERT_NE(&arr[i], &copy_arr[i]);
  });

  arr.free();
  copy_arr.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, DeepCopyPolygon) {
  using WORK_EXEC_POLICY = TypeParam;

  using MV = Spheral::ManagedVector<Spheral::GeomPolygon>;
  MV arr(6);

  RAJA::forall<WORK_EXEC_POLICY>(
      TRS_UINT(0, arr.size()),
      [=] RAJA_HOST_DEVICE(unsigned i) 
    { arr[i] = Spheral::GeomPolygon({{0,0},{1,0},{1,1},{0.5,1.5},{0,1}}); });

  SPHERAL_ASSERT_EQ(arr.size(), 6);
  MV copy_arr = arr.clone();

  SPHERAL_ASSERT_NE(&arr[0], &copy_arr[0]);

  SPHERAL_ASSERT_EQ(arr.size(), copy_arr.size());
  SPHERAL_ASSERT_EQ(arr.capacity(), copy_arr.capacity());

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 6), [=] RAJA_HOST(unsigned i) {
    SPHERAL_ASSERT_EQ(copy_arr[i], arr[i]);
    SPHERAL_ASSERT_NE(&arr[i], &copy_arr[i]);
  });

  copy_arr.free();
  arr.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, AssignmentOperator) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr(6, 5);

  MVDouble copy_arr = arr;

  RAJA::forall<WORK_EXEC_POLICY>(
      TRS_UINT(0, 6), [=] RAJA_HOST_DEVICE(unsigned i) { arr[i] = i; });

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 6), [=] RAJA_HOST(unsigned i) {
    SPHERAL_ASSERT_EQ(copy_arr[i], i);
  });

  SPHERAL_ASSERT_EQ(&arr[0], &copy_arr[0]);

  copy_arr.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, Equivalence) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr(6, 5);
  MVDouble arr2(6, 3);

  MVDouble copy_arr = arr;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)

  for (size_t i = 0; i < 6u; i++) {
    arr[i] = i;
    arr2[i] = i;
  }
  SPHERAL_ASSERT_TRUE(copy_arr == arr);
  SPHERAL_ASSERT_FALSE(copy_arr == arr2);
  EXEC_IN_SPACE_END()

  EXEC_IN_SPACE_BEGIN(LOOP_EXEC_POLICY)
  ASSERT_TRUE(copy_arr == arr);
  ASSERT_FALSE(copy_arr == arr2);
  EXEC_IN_SPACE_END()

  arr2.free();
  copy_arr.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, PushBackDefault) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr;
  double val = 5;

  arr.push_back(val);

  SPHERAL_ASSERT_EQ(arr.size(), 1u);
  SPHERAL_ASSERT_EQ(arr.capacity(), MVDouble::initial_capacity);
  SPHERAL_ASSERT_EQ(arr[0], val);

  MVDouble arr2;

  for (size_t i = 0; i < 10; i++)
    arr2.push_back(i);

  for (size_t i = 0; i < arr2.size(); i++)
    SPHERAL_ASSERT_EQ(arr2[i], i);

  RAJA::forall<WORK_EXEC_POLICY>(TRS_UINT(0, 10),
                                 [=] RAJA_HOST_DEVICE(unsigned idx) {
                                   SPHERAL_ASSERT_EQ(arr2.size(), 10u);
                                   SPHERAL_ASSERT_EQ(arr2[idx], idx);
                                   arr2[idx] *= 2;
                                 });

  RAJA::forall<LOOP_EXEC_POLICY>(TRS_UINT(0, 6), [=] RAJA_HOST(unsigned idx) {
    SPHERAL_ASSERT_EQ(arr2[idx], idx * 2);
  });

  SPHERAL_ASSERT_EQ(arr2.capacity(), MVDouble::initial_capacity * 2);

  arr.free();
  arr2.free();
}

TEST(ManagedVectorTest, PushBackMove) {
  MVDouble arr;

  arr.push_back(std::move(5));

  SPHERAL_ASSERT_EQ(arr.size(), 1u);
  SPHERAL_ASSERT_EQ(arr.capacity(), MVDouble::initial_capacity);
  SPHERAL_ASSERT_EQ(arr[0], 5);

  MVDouble arr2;

  for (size_t i = 0; i < 10; i++)
    arr2.push_back(std::move(i));

  for (size_t i = 0; i < arr2.size(); i++)
    SPHERAL_ASSERT_EQ(arr2[i], i);

  SPHERAL_ASSERT_EQ(arr2.capacity(), MVDouble::initial_capacity * 2);

  arr.free();
  arr2.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, ResizeLargerNoRealloc) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr.size(), 0);
  SPHERAL_ASSERT_EQ(arr.capacity(), 0);
  EXEC_IN_SPACE_END()

  arr.move(chai::CPU);
  arr.resize(10);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr.size(), 10);
  SPHERAL_ASSERT_EQ(arr.capacity(), 16);
  EXEC_IN_SPACE_END()

  MVDouble arr2(4);
  MVValType *pre_resize_ptr = &(arr2[0]);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr2.size(), 4);
  // SPHERAL_ASSERT_EQ(arr2.capacity(), MVDouble::initial_capacity);
  EXEC_IN_SPACE_END()

  arr.move(chai::CPU);
  arr2.resize(6);
  MVValType *post_resize_ptr = &(arr2[0]);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr2.size(), 6);
  // SPHERAL_ASSERT_EQ(arr2.capacity(), MVDouble::initial_capacity);
  EXEC_IN_SPACE_END()

  // No reallocation has been performed the data should be in the same location.
  ASSERT_EQ(pre_resize_ptr, post_resize_ptr);

  arr.free();
  arr2.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, ResizeLargerRealloc) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr(4);
  MVValType *pre_resize_ptr = &(arr[0]);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr.size(), 4);
  // SPHERAL_ASSERT_EQ(arr.capacity(), 0);
  EXEC_IN_SPACE_END()

  arr.move(chai::CPU);
  arr.resize(12);
  MVValType *post_resize_ptr = &(arr[0]);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr.size(), 12);
  SPHERAL_ASSERT_EQ(arr.capacity(), 16);
  EXEC_IN_SPACE_END()

  // Reallocation has been performed the data should NOT be in the same
  // location.
  ASSERT_NE(pre_resize_ptr, post_resize_ptr);

  arr.free();
}

GPU_TYPED_TEST(ManagedVectorTypedTest, ResizeSmaller) {
  using WORK_EXEC_POLICY = TypeParam;

  MVDouble arr(4);
  MVValType *pre_resize_ptr = &(arr[0]);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr.size(), 4);
  SPHERAL_ASSERT_EQ(arr.capacity(), MVDouble::initial_capacity);
  EXEC_IN_SPACE_END()

  arr.resize(2);
  MVValType *post_resize_ptr = &(arr[0]);

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  SPHERAL_ASSERT_EQ(arr.size(), 2);
  SPHERAL_ASSERT_EQ(arr.capacity(), MVDouble::initial_capacity);
  EXEC_IN_SPACE_END()

  // No reallocation has been performed the data should be in the same location.
  ASSERT_EQ(pre_resize_ptr, post_resize_ptr);

  arr.free();
}

TEST(ManagedVectorTest, Erase) {
  MVDouble arr;
  std::vector<double> check = {0, 1, 2, 3, 4, 5};

  for (size_t i = 0; i < 6; i++) {
    arr.push_back(check[i]);
  }
  SPHERAL_ASSERT_EQ(arr.size(), check.size());

  // Erase the last element
  arr.erase(arr.end() - 1);

  std::vector<double> check2 = {0, 1, 2, 3, 4};
  SPHERAL_ASSERT_EQ(arr.size(), check2.size());
  for (size_t i = 0; i < arr.size(); i++)
    SPHERAL_ASSERT_EQ(arr[i], check2[i]);

  // Erase the 3rd element
  arr.erase(arr.begin() + 2);

  std::vector<double> check3 = {0, 1, 3, 4};
  SPHERAL_ASSERT_EQ(arr.size(), check3.size());
  for (size_t i = 0; i < arr.size(); i++)
    SPHERAL_ASSERT_EQ(arr[i], check3[i]);

  // Erase the first element
  arr.erase(arr.begin());

  std::vector<double> check4 = {1, 3, 4};
  SPHERAL_ASSERT_EQ(arr.size(), check4.size());
  for (size_t i = 0; i < arr.size(); i++)
    SPHERAL_ASSERT_EQ(arr[i], check4[i]);

  arr.free();
}

TEST(ManagedVectorTest, Insert) {
  MVDouble arr;
  std::vector<double> check = {0, 1, 2, 3, 4, 5};

  for (size_t i = 0; i < 6; i++) {
    arr.insert(arr.begin() + i, check[i]);
  }
  SPHERAL_ASSERT_EQ(arr.size(), check.size());
  for (size_t i = 0; i < arr.size(); i++)
    SPHERAL_ASSERT_EQ(arr[i], check[i]);

  // Insert element at end
  arr.insert(arr.end(), 6);

  std::vector<double> check2 = {0, 1, 2, 3, 4, 5, 6};
  SPHERAL_ASSERT_EQ(arr.size(), check2.size());
  for (size_t i = 0; i < arr.size(); i++)
    SPHERAL_ASSERT_EQ(arr[i], check2[i]);

  // Erase the 3rd element
  arr.insert(arr.begin() + 2, 7);

  std::vector<double> check3 = {0, 1, 7, 2, 3, 4, 5, 6};
  SPHERAL_ASSERT_EQ(arr.size(), check3.size());
  for (size_t i = 0; i < arr.size(); i++)
    SPHERAL_ASSERT_EQ(arr[i], check3[i]);

  // Erase the first element
  arr.insert(arr.begin(), -1);

  std::vector<double> check4 = {-1, 0, 1, 7, 2, 3, 4, 5, 6};
  SPHERAL_ASSERT_EQ(arr.size(), check4.size());
  for (size_t i = 0; i < arr.size(); i++)
    SPHERAL_ASSERT_EQ(arr[i], check4[i]);

  arr.free();
}

TEST(ManagedVectorTest, DeepCopy) {
  double init_val = 5;
  MVDouble arr(6, init_val);

  MVDouble copy_arr(arr);
  MVDouble deep_copy_arr(deepCopy(arr));

  for (size_t i = 0; i < 6u; i++) {
    arr[i] = 4;
  }

  for (size_t i = 0; i < 6u; i++) {
    ASSERT_NE(arr[i], deep_copy_arr[i]);
    SPHERAL_ASSERT_EQ(deep_copy_arr[i], init_val);
  }

  ASSERT_FALSE(arr == deep_copy_arr);
  ASSERT_NE(&arr[0], &deep_copy_arr[0]);

  copy_arr.free();
  deep_copy_arr.free();
}
