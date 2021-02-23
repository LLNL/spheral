#include "Field/Field.hh"
#include "Utilities/SidreDataCollection.hh"
#include "Geometry/Dimension.hh"
#include "axom/sidre.hpp"

#include "gtest/gtest.h"
#include <iostream>

// #include <typeinfo>
// #include <vector>

template <typename T>
class SidreDataCollectionTest : public ::testing::Test
{
  public:
    Spheral::SidreDataCollection myData;
    int n = 100;
    T* rawSidreData;

    Spheral::Field<Spheral::Dim<1>, T> makeField()
    {
      Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
      Spheral::Field<Spheral::Dim<1>, T> testField("test field", makeNodeList);
      for (int i = 0; i < n; i++)
        testField[i] = i;
      return testField;
    }

    void allocRawSidreData(const Spheral::Field<Spheral::Dim<1>, T>& testField)
    {
      rawSidreData = myData.alloc_view("SidreTest", testField)->getData();
    }
};

using MyTypes = ::testing::Types<char, int, size_t, uint32_t, uint64_t, float, double>;
TYPED_TEST_SUITE(SidreDataCollectionTest, MyTypes);

TYPED_TEST(SidreDataCollectionTest, scalar)
{
  auto testField = this->makeField();
  
  // for (int i = 0; i < this->n; i++)
  //   std::cout << testField[i] << " ";
  // std::cout << "\n";

  this->allocRawSidreData(testField);
  
  //this->myData.printDataStore();

  // std::cout << "Result: ";
  // for (int i = 0; i < this->n; i++)
  //   std::cout << this->rawSidreData[i] << " ";
  // std::cout << "\n";
  
  for (int i = 0; i < this->n; i++)
    EXPECT_EQ(testField[i], this->rawSidreData[i]);
}



TEST(SidreDataCollectionTestString, string)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, std::string> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = "This is a test string: " + std::to_string(i) + "\n";

  // for (int i = 0; i < n; i++)
  //   std::cout << testField[i];
  
  char *rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  // myData.printDataStore();

  for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0][i], rawSidreData[i]);
}



template <typename T>
class SidreDataCollectionTestVector : public ::testing::Test
{
  public:
    Spheral::SidreDataCollection myData;
    int n = 10;
    T* rawSidreData;

    Spheral::Field<Spheral::Dim<1>, std::vector<T>> makeField()
    {
      Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
      Spheral::Field<Spheral::Dim<1>, std::vector<T>> testField("test field", makeNodeList);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          testField[i].emplace_back(i);
      return testField;
    }

    void allocRawSidreData(const Spheral::Field<Spheral::Dim<1>, std::vector<T>>& testField)
    {
      rawSidreData = myData.alloc_view("SidreTest", testField)->getData();
    }
};

TYPED_TEST_SUITE(SidreDataCollectionTestVector, MyTypes);

TYPED_TEST(SidreDataCollectionTestVector, vector)
{
  auto testField = this->makeField();

  // for (int i = 0; i < n; i++)
  // {
  //   for (int j = 0; j < n; j++)
  //     std::cout << testField[i][j] << " ";
  //   std::cout << std::endl;
  // }

  this->allocRawSidreData(testField);

  //this->myData.printDataStore();

  for (int i = 0; i < this->n; i++)
    EXPECT_EQ(testField[0][i], this->rawSidreData[i]);
}



template <typename T>
class SidreDataCollectionTestTupleThree : public ::testing::Test
{
  public:
    Spheral::SidreDataCollection myData;
    int n = 5;
    T* rawSidreData;

    Spheral::Field<Spheral::Dim<1>, std::tuple<T, T, T>> makeField()
    {
      Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
      Spheral::Field<Spheral::Dim<1>, std::tuple<T, T, T>> testField("test field", makeNodeList);
      for (int i = 0; i < n; i++)
        testField[i] = std::make_tuple(i, i, i);
      return testField;
    }

    void allocRawSidreData(const Spheral::Field<Spheral::Dim<1>, std::tuple<T, T, T>>& testField)
    {
      rawSidreData = myData.alloc_view("SidreTest", testField)->getData();
    }
};

TYPED_TEST_SUITE(SidreDataCollectionTestTupleThree, MyTypes);

TYPED_TEST(SidreDataCollectionTestTupleThree, tuple)
{
  auto testField = this->makeField();

  // for (int i = 0; i < this->n; i++)
  //   std::cout << testField[i] << " ";
  // std::cout << std::endl;

  this->allocRawSidreData(testField);

  // this->myData.printDataStore();

  EXPECT_EQ(std::get<0>(testField[0]), this->rawSidreData[0]);
  EXPECT_EQ(std::get<1>(testField[0]), this->rawSidreData[1]);
  EXPECT_EQ(std::get<2>(testField[0]), this->rawSidreData[2]);
}

int main(int argc, char** argv)
{
// #ifdef AXOM_USE_MPI
//   MPI_Init(&argc, &argv);
// #endif
  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

// #ifdef AXOM_USE_MPI
//   MPI_Finalize();
// #endif

  return result;
}

