#include "Field/Field.hh"
#include "Utilities/SidreDataCollection.hh"
#include "Geometry/Dimension.hh"
#include "axom/sidre.hpp"

#include "gtest/gtest.h"
#include <iostream>

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

  this->allocRawSidreData(testField);
  
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
  
  char *rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

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

  this->allocRawSidreData(testField);

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

  this->allocRawSidreData(testField);

  EXPECT_EQ(std::get<0>(testField[0]), this->rawSidreData[0]);
  EXPECT_EQ(std::get<1>(testField[0]), this->rawSidreData[1]);
  EXPECT_EQ(std::get<2>(testField[0]), this->rawSidreData[2]);
}



TEST(SidreDataCollectionTestDim1Vector, Dim1Vector)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<1>::Vector(i);
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  myData.printDataStore();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0].x(), rawSidreData[1]);
}

TEST(SidreDataCollectionTestDim1Vector3d, Dim1Vector3d)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector3d> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<1>::Vector3d(i, i + 1, i + 2);

  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //for (int i = 0; i < n; i++)
  {
    EXPECT_EQ(testField[0].x(), rawSidreData[0]);
    EXPECT_EQ(testField[0].y(), rawSidreData[1]);
    EXPECT_EQ(testField[0].z(), rawSidreData[2]);
  }
}

TEST(SidreDataCollectionTestDim1Tensor, Tensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Tensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<1>::Tensor(i);

  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0].xx(), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim1SymTensor, SymTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<1>::SymTensor(i);

  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0].xx(), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim1ThirdRankTensor, ThirdRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::ThirdRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<1>::ThirdRankTensor(i);

  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0,0,0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim1FourthRankTensor, FourthRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::FourthRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<1>::FourthRankTensor(i);

  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0, 0, 0, 0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim1FifthRankTensor, FifthRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<1>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::FifthRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<1>::FifthRankTensor(i);

  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0, 0, 0, 0, 0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim2Vector, Dim2Vector)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<2>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Vector> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<2>::Vector(i, i + 1);

  // std::cout << testField[0].x()
  //           << testField[0].y() << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
  {
    EXPECT_EQ(testField[0].x(), rawSidreData[0]);
    EXPECT_EQ(testField[0].y(), rawSidreData[1]);
  }
}

TEST(SidreDataCollectionTestDim2Tensor, Dim2Tensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<2>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Tensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<2>::Tensor(i, i + 1, i + 2, i + 3);

  // std::cout << testField[0].xx()
  //           << testField[0].xy() << std::endl
  //           << testField[0].yx()
  //           << testField[0].yy() << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
  {
    EXPECT_EQ(testField[0].xx(), rawSidreData[0]);
    EXPECT_EQ(testField[0].xy(), rawSidreData[1]);
    EXPECT_EQ(testField[0].yx(), rawSidreData[2]);
    EXPECT_EQ(testField[0].yy(), rawSidreData[3]);
  }
}

TEST(SidreDataCollectionTestDim2SymTensor, SymTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<2>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<2>::SymTensor(i, i + 1, i + 1, i);

  // std::cout << testField[0].xx()
  //           << testField[0].xy() << std::endl
  //           << testField[0].yx()
  //           << testField[0].yy() << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  EXPECT_EQ(testField[0].xx(), rawSidreData[0]);
  EXPECT_EQ(testField[0].xy(), rawSidreData[1]);
}

TEST(SidreDataCollectionTestDim2ThirdRankTensor, ThirdRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<2>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::ThirdRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<2>::ThirdRankTensor(i);

  //std::cout << testField[0](0, 0, 0) << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0,0,0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim2FourthRankTensor, FourthRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<2>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::FourthRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<2>::FourthRankTensor(i);

  //std::cout << testField[0](0, 0, 0, 0) << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0, 0, 0, 0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim2FifthRankTensor, FifthRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<2>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::FifthRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<2>::FifthRankTensor(i);

  //std::cout << testField[0](0, 0, 0, 0, 0) << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0, 0, 0, 0, 0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim3Vector, Dim3Vector)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<3>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Vector> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<3>::Vector(i, i + 1, i + 2);

  // std::cout << testField[0].x()
  //           << testField[0].y() << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
  {
    EXPECT_EQ(testField[0].x(), rawSidreData[0]);
    EXPECT_EQ(testField[0].y(), rawSidreData[1]);
    EXPECT_EQ(testField[0].z(), rawSidreData[2]);
  }
}

TEST(SidreDataCollectionTestDim3Tensor, Dim3Tensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<3>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Tensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<3>::Tensor(i, i + 1, i + 2,
                                           i, i + 1, i + 2,
                                           i, i + 1, i + 2);

  // std::cout << testField[0].xx()
  //           << testField[0].xy() << std::endl
  //           << testField[0].yx()
  //           << testField[0].yy() << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
  {
    EXPECT_EQ(testField[0].xx(), rawSidreData[0]);
    EXPECT_EQ(testField[0].xy(), rawSidreData[1]);
    EXPECT_EQ(testField[0].xz(), rawSidreData[2]);
    EXPECT_EQ(testField[0].yx(), rawSidreData[3]);
    EXPECT_EQ(testField[0].yy(), rawSidreData[4]);
    EXPECT_EQ(testField[0].yz(), rawSidreData[5]);
    EXPECT_EQ(testField[0].zx(), rawSidreData[6]);
    EXPECT_EQ(testField[0].zy(), rawSidreData[7]);
    EXPECT_EQ(testField[0].zz(), rawSidreData[8]);
  }
}

TEST(SidreDataCollectionTestDim3SymTensor, SymTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<3>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<3>::SymTensor(i, i + 1, i + 2,
                                              i + 1, i, i + 1,
                                              i + 2, i + 1, i);

  // std::cout << testField[0].xx()
  //           << testField[0].xy()
  //           << testField[0].xz() << std::endl
  //           << testField[0].yx()
  //           << testField[0].yy()
  //           << testField[0].yz() << std::endl
  //           << testField[0].zx()
  //           << testField[0].zy()
  //           << testField[0].zz() << std::endl;

  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  EXPECT_EQ(testField[0].xx(), rawSidreData[0]);
  EXPECT_EQ(testField[0].xy(), rawSidreData[1]);
  EXPECT_EQ(testField[0].xz(), rawSidreData[2]);
}

TEST(SidreDataCollectionTestDim3ThirdRankTensor, ThirdRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<3>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::ThirdRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<3>::ThirdRankTensor(i);

  //std::cout << testField[0](0, 0, 0) << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0,0,0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim3FourthRankTensor, FourthRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<3>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::FourthRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<3>::FourthRankTensor(i);

  //std::cout << testField[0](0, 0, 0, 0) << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0, 0, 0, 0), rawSidreData[0]);
}

TEST(SidreDataCollectionTestDim3FifthRankTensor, FifthRankTensor)
{
  Spheral::SidreDataCollection myData;
  int n = 10;

  Spheral::NodeList<Spheral::Dim<3>> makeNodeList("test bed", n, 0);
  Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::FifthRankTensor> testField("test field", makeNodeList);
  for (int i = 0; i < n; i++)
    testField[i] = Spheral::Dim<3>::FifthRankTensor(i);

  //std::cout << testField[0](0, 0, 0, 0, 0) << std::endl;
  
  double* rawSidreData = myData.alloc_view("SidreTest", testField)->getData();

  //myData.printDataStore();

  //for (int i = 0; i < n; i++)
    EXPECT_EQ(testField[0](0, 0, 0, 0, 0), rawSidreData[0]);
}


// int main(int argc, char** argv)
// {
// // #ifdef AXOM_USE_MPI
// //   MPI_Init(&argc, &argv);
// // #endif
//   ::testing::InitGoogleTest(&argc, argv);

//   int result = RUN_ALL_TESTS();

// // #ifdef AXOM_USE_MPI
// //   MPI_Finalize();
// // #endif

//   return result;
// }

