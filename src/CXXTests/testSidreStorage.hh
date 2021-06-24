//---------------------------------Spheral++----------------------------------//
// testSidreStorage -- Unit tests for SidreDataCollection
//
//
// Created by Mikhail Zakharchanka, 2020
//----------------------------------------------------------------------------//

#include "Field/Field.hh"
#include "Utilities/SidreDataCollection.hh"
#include "Geometry/Dimension.hh"
#include "axom/sidre.hpp"

#include "gtest/gtest.h"

using arithmeticTypes = ::testing::Types<char, int, size_t, uint32_t, uint64_t, float, double>;

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
      for (int i = 0; i < n; ++i)
        testField[i] = i;
      return testField;
    }

    void allocRawSidreData(const Spheral::Field<Spheral::Dim<1>, T>& testField)
    {
      rawSidreData = myData.alloc_view("SidreTest", testField)->getData();
    }
};

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
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; j++)
          testField[i].emplace_back(i);
      return testField;
    }

    void allocRawSidreData(const Spheral::Field<Spheral::Dim<1>, std::vector<T>>& testField)
    {
      rawSidreData = myData.alloc_view("SidreTest", testField)->getData();
    }
};

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
      for (int i = 0; i < n; ++i)
        testField[i] = std::make_tuple(i, i, i);
      return testField;
    }

    void allocRawSidreData(const Spheral::Field<Spheral::Dim<1>, std::tuple<T, T, T>>& testField)
    {
      rawSidreData = myData.alloc_view("SidreTest", testField)->getData();
    }
};


TYPED_TEST_SUITE(SidreDataCollectionTestTupleThree, arithmeticTypes);
TYPED_TEST_SUITE(SidreDataCollectionTest, arithmeticTypes);
TYPED_TEST_SUITE(SidreDataCollectionTestVector, arithmeticTypes);
