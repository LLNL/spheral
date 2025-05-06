#include "Utilities/ManagedVector.hh"
#include "test-utilities.hh"

// --- Test Data Structure (same as before) ---
struct TestPayload {
  int id;
  std::string data;
  double value;

  TestPayload(int i = 0, std::string d = "default", double v = 0.0)
      : id(i), data(std::move(d)), value(v) {}

  TestPayload(const TestPayload &other)
      : id(other.id), data(other.data), value(other.value) {}

  TestPayload(TestPayload &&other) noexcept
      : id(other.id), data(std::move(other.data)), value(other.value) {
    other.id = -1;
  }

  TestPayload &operator=(const TestPayload &other) {
    if (this != &other) {
      id = other.id;
      data = other.data;
      value = other.value;
    }
    return *this;
  }

  TestPayload &operator=(TestPayload &&other) noexcept {
    if (this != &other) {
      id = other.id;
      data = std::move(other.data);
      value = other.value;
      other.id = -1;
    }
    return *this;
  }

  ~TestPayload() =
      default; // Destructor output removed for cleaner benchmark output

  bool operator==(const TestPayload &other) const {
    return id == other.id && data == other.data && value == other.value;
  }
};

// --- Timer Utility (same as before) ---
class Timer {
public:
  Timer(const std::string &name)
      : name_(name), start_time_(std::chrono::high_resolution_clock::now()) {}
  ~Timer() { stop(); }
  void stop() {
    if (!stopped_) {
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
          end_time - start_time_);
      std::cout << "[TIMER] " << name_ << " took " << duration.count() << " us."
                << std::endl;
      stopped_ = true;
    }
  }

private:
  std::string name_;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
  bool stopped_ = false;
};

// --- Global constants for test sizes ---
const size_t NUM_ELEMENTS_SMALL = 10000;
const size_t NUM_ELEMENTS_MEDIUM =
    100000; // Reduced for faster gtest runs, adjust as needed
const size_t NUM_ELEMENTS_LARGE = 1000000; // Reduced for faster gtest runs

// --- Test Fixture for Typed Tests ---
template <typename VectorTypeParam>
class VectorPerformanceTest : public ::testing::Test {
public:
  using VecType = VectorTypeParam; // TypeParam is the vector type (e.g.
                                   // MyVector<TestPayload>)

  std::string getVectorName() {
    // Note: MyVector<TestPayload, std::allocator<TestPayload>> is the full type
    // This check is simplified. For more complex types, you might need more
    // robust type identification.
    if constexpr (std::is_same_v<VecType,
                                 Spheral::ManagedVector<TestPayload>> ||
                  std::is_same_v<
                      VecType,
                      Spheral::ManagedVector<TestPayload>>) { // Allow
                                                              // default
                                                              // allocator
      return "Spheral::ManagedVector";

    } else if constexpr (std::is_same_v<
                             VecType,
                             std::vector<TestPayload,
                                         std::allocator<TestPayload>>> ||
                         std::is_same_v<VecType, std::vector<TestPayload>>) {
      return "std::vector";
    }

    return "UnknownVectorType";
  }
};

// Define the types to be tested
// Using MyVector<TestPayload> assumes it uses std::allocator<TestPayload> by
// default If your MyVector requires explicit allocator, use
// MyVector<TestPayload, std::allocator<TestPayload>>
using VectorTypes = ::testing::Types<
    Spheral::ManagedVector<TestPayload>, // Or MyVector<TestPayload,
                                         // std::allocator<TestPayload>>
    std::vector<TestPayload>>;
TYPED_TEST_SUITE(VectorPerformanceTest, VectorTypes);

// --- Typed Test Cases ---

TYPED_TEST(VectorPerformanceTest, PushBackLValue) {
  using VecType = typename TestFixture::VecType;
  VecType vec;
  std::string vec_name = TestFixture::getVectorName();
  Timer t(vec_name + " PushBackLValue (" + std::to_string(NUM_ELEMENTS_MEDIUM) +
          " elements)");
  for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
    vec.push_back(TestPayload(static_cast<int>(i),
                              "data_item_" + std::to_string(i),
                              static_cast<double>(i)));
  }
  EXPECT_EQ(vec.size(), NUM_ELEMENTS_MEDIUM);
}

TYPED_TEST(VectorPerformanceTest, PushBackRValue) {
  using VecType = typename TestFixture::VecType;
  VecType vec;
  std::string vec_name = TestFixture::getVectorName();
  Timer t(vec_name + " PushBackRValue (" + std::to_string(NUM_ELEMENTS_MEDIUM) +
          " elements)");
  for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
    vec.push_back(TestPayload(static_cast<int>(i),
                              "move_data_item_" + std::to_string(i),
                              static_cast<double>(i) + 0.5));
  }
  EXPECT_EQ(vec.size(), NUM_ELEMENTS_MEDIUM);
}

TYPED_TEST(VectorPerformanceTest, EmplaceBack) {
  using VecType = typename TestFixture::VecType;
  VecType vec;
  std::string vec_name = TestFixture::getVectorName();
  Timer t(vec_name + " EmplaceBack (" + std::to_string(NUM_ELEMENTS_MEDIUM) +
          " elements)");
  for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
    vec.emplace_back(static_cast<int>(i), "emplace_data_" + std::to_string(i),
                     static_cast<double>(i) + 0.75);
  }
  EXPECT_EQ(vec.size(), NUM_ELEMENTS_MEDIUM);
}

TYPED_TEST(VectorPerformanceTest, ConstructNWithValue) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  TestPayload default_val(0, "construct_val", 1.0);
  Timer t(vec_name + " ConstructNWithValue (" +
          std::to_string(NUM_ELEMENTS_MEDIUM) + " elements)");
  VecType vec(NUM_ELEMENTS_MEDIUM, default_val);
  EXPECT_EQ(vec.size(), NUM_ELEMENTS_MEDIUM);
  // if (!vec.empty()) {
  EXPECT_EQ(vec[0], default_val); // Basic check
  //}
}

// Note: MyVector stub doesn't have default constructor with count. std::vector
// does. Add this test if your MyVector supports `MyVector(count)` for default
// construction.

TYPED_TEST(VectorPerformanceTest, ConstructNDefault) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  if (vec_name == "MyVector" && !std::is_constructible_v<VecType, size_t>) {
    GTEST_SKIP() << "MyVector(count) for default construction not implemented "
                 << "in stub or full class.";
  }
  Timer t(vec_name + " ConstructNDefault (" +
          std::to_string(NUM_ELEMENTS_MEDIUM) + " elements)");
  VecType vec(
      NUM_ELEMENTS_MEDIUM); // Assumes default construction of TestPayload
  EXPECT_EQ(vec.size(), NUM_ELEMENTS_MEDIUM);
}

/*
TYPED_TEST(VectorPerformanceTest, CopyConstructor) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  VecType original_vec;
  for (size_t i = 0; i < NUM_ELEMENTS_SMALL; ++i) { // Smaller for copy
    original_vec.emplace_back(static_cast<int>(i),
                              "original_item_" + std::to_string(i),
                              static_cast<double>(i));
  }
  Timer t(vec_name + " CopyConstructor (" + std::to_string(NUM_ELEMENTS_SMALL) +
          " elements)");
  VecType copied_vec = original_vec;
  EXPECT_EQ(copied_vec.size(), original_vec.size());
  if (!original_vec.empty() && !copied_vec.empty()) {
    EXPECT_EQ(copied_vec[0], original_vec[0]);
  }
}

TYPED_TEST(VectorPerformanceTest, CopyAssignment) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  VecType original_vec;
  for (size_t i = 0; i < NUM_ELEMENTS_SMALL; ++i) { // Smaller for copy
    original_vec.emplace_back(static_cast<int>(i),
                              "original_item_" + std::to_string(i),
                              static_cast<double>(i));
  }
  VecType assigned_vec;
  Timer t(vec_name + " CopyAssignment (" + std::to_string(NUM_ELEMENTS_SMALL) +
          " elements)");
  assigned_vec = original_vec;
  EXPECT_EQ(assigned_vec.size(), original_vec.size());
  if (!original_vec.empty() && !assigned_vec.empty()) {
    EXPECT_EQ(assigned_vec[0], original_vec[0]);
  }
}

TYPED_TEST(VectorPerformanceTest, MoveConstructor) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  auto create_vec = [&]() {
    VecType vec;
    for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
      vec.emplace_back(static_cast<int>(i),
                       "move_original_item_" + std::to_string(i),
                       static_cast<double>(i));
    }
    return vec;
  };
  VecType source_vec = create_vec();
  size_t original_size = source_vec.size();
  Timer t(vec_name + " MoveConstructor (" +
          std::to_string(NUM_ELEMENTS_MEDIUM) + " elements)");
  VecType moved_vec = std::move(source_vec);
  EXPECT_EQ(moved_vec.size(), original_size);
  // Per standard, moved-from std::vector is in a valid but unspecified state.
  // For MyVector, we expect it to be empty if implemented correctly.
  if (vec_name == "MyVector") {
    EXPECT_TRUE(source_vec.empty()); // Or check size == 0
  }
}

TYPED_TEST(VectorPerformanceTest, MoveAssignment) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  auto create_vec = [&]() {
    VecType vec;
    for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
      vec.emplace_back(static_cast<int>(i),
                       "move_original_item_" + std::to_string(i),
                       static_cast<double>(i));
    }
    return vec;
  };
  VecType source_vec = create_vec();
  size_t original_size = source_vec.size();
  VecType moved_assign_vec;
  Timer t(vec_name + " MoveAssignment (" + std::to_string(NUM_ELEMENTS_MEDIUM) +
          " elements)");
  moved_assign_vec = std::move(source_vec);
  EXPECT_EQ(moved_assign_vec.size(), original_size);
  if (vec_name == "MyVector") {
    EXPECT_TRUE(source_vec.empty());
  }
}
*/

TYPED_TEST(VectorPerformanceTest, ReserveThenPushBack) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  VecType vec;
  {
    Timer t(vec_name + " Reserve (" + std::to_string(NUM_ELEMENTS_MEDIUM) +
            " elements)");
    vec.reserve(NUM_ELEMENTS_MEDIUM);
  }
  EXPECT_GE(vec.capacity(), NUM_ELEMENTS_MEDIUM);
  {
    Timer t(vec_name + " PushBackAfterReserve (" +
            std::to_string(NUM_ELEMENTS_MEDIUM) + " elements)");
    for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
      vec.push_back(TestPayload(static_cast<int>(i),
                                "reserved_item_" + std::to_string(i),
                                static_cast<double>(i)));
    }
  }
  EXPECT_EQ(vec.size(), NUM_ELEMENTS_MEDIUM);
}

TYPED_TEST(VectorPerformanceTest, AccessOperatorBracket) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  if (NUM_ELEMENTS_SMALL == 0)
    GTEST_SKIP() << "Skipping access test for 0 elements.";
  VecType vec;
  for (size_t i = 0; i < NUM_ELEMENTS_SMALL; ++i) {
    vec.emplace_back(static_cast<int>(i));
  }
  volatile int sink = 0;
  Timer t(vec_name + " Operator[] Access (" +
          std::to_string(NUM_ELEMENTS_LARGE) + " accesses)");
  for (size_t i = 0; i < NUM_ELEMENTS_LARGE; ++i) {
    sink += vec[i % NUM_ELEMENTS_SMALL].id;
  }
  (void)sink; // Use sink
}

/*
TYPED_TEST(VectorPerformanceTest, AccessAtMethod) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  if (NUM_ELEMENTS_SMALL == 0)
    GTEST_SKIP() << "Skipping access test for 0 elements.";
  VecType vec;
  for (size_t i = 0; i < NUM_ELEMENTS_SMALL; ++i) {
    vec.emplace_back(static_cast<int>(i));
  }
  volatile int sink = 0;
  Timer t(vec_name + " At() Access (" + std::to_string(NUM_ELEMENTS_LARGE) +
          " accesses)");
  for (size_t i = 0; i < NUM_ELEMENTS_LARGE; ++i) {
    sink += vec.at(i % NUM_ELEMENTS_SMALL).id;
  }
  (void)sink;
}
*/

/*
TYPED_TEST(VectorPerformanceTest, PopBack) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  VecType vec;
  for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
    vec.emplace_back(static_cast<int>(i));
  }
  Timer t(vec_name + " PopBack (" + std::to_string(NUM_ELEMENTS_MEDIUM) +
          " elements)");
  for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
    vec.pop_back();
  }
  EXPECT_TRUE(vec.empty());
}
*/

TYPED_TEST(VectorPerformanceTest, Clear) {
  using VecType = typename TestFixture::VecType;
  std::string vec_name = TestFixture::getVectorName();
  VecType vec;
  for (size_t i = 0; i < NUM_ELEMENTS_MEDIUM; ++i) {
    vec.emplace_back(static_cast<int>(i));
  }
  Timer t(vec_name + " Clear (" + std::to_string(NUM_ELEMENTS_MEDIUM) +
          " elements)");
  vec.clear();
  // EXPECT_TRUE(vec.empty());
}

// --- Main function to run gtests ---
int main(int argc, char **argv) {
  std::cout << "Starting MyVector vs std::vector Performance Comparison (gtest)"
            << std::endl;
  std::cout << "=============================================================="
            << std::endl;
  std::cout << "Element Type: TestPayload (int, std::string, double)"
            << std::endl;
  std::cout << "Compiling with optimizations (e.g., -O2 or -O3) is crucial for "
               "meaningful results."
            << std::endl;
  std::cout << "Ensure MyVector.hpp is correctly included and linked."
            << std::endl;
  std::cout << "Number of elements for most tests (MEDIUM): "
            << NUM_ELEMENTS_MEDIUM << std::endl;
  std::cout << "Number of elements for copy tests (SMALL): "
            << NUM_ELEMENTS_SMALL << std::endl;
  std::cout << "Number of accesses for access tests (LARGE): "
            << NUM_ELEMENTS_LARGE << std::endl;

  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  std::cout
      << "\n=============================================================="
      << std::endl;
  std::cout << "Benchmark finished. Check gtest output and [TIMER] lines for "
               "performance."
            << std::endl;
  return result;
}
