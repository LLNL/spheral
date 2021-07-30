//---------------------------------Spheral++----------------------------------//
// testSidreStorage -- Unit tests for SidreDataCollection
//
//
// Created by Mikhail Zakharchanka, 2021
//----------------------------------------------------------------------------//

#include "Field/Field.hh"
#include "Utilities/SidreDataCollection.hh"
#include "Geometry/Dimension.hh"
#include "axom/sidre.hpp"
#include "gtest/gtest.h"


template<int Dim, typename T>
struct SpheralTypeInfo {
  using DIM_T = Spheral::Dim<Dim>;
  using SPHERAL_T = T;
}; 

template<typename T_DATA>
using SpheralTestField = Spheral::Field<typename T_DATA::DIM_T, typename T_DATA::SPHERAL_T> ;
template<typename T_DATA>
using SpheralTestNodeList = Spheral::NodeList<typename T_DATA::DIM_T> ;

#define INIT_FIELD( _dim, _type)          \
template<typename T_DATA,                 \
         typename std::enable_if<         \
           std::is_same<                  \
             T_DATA,                      \
             SpheralTypeInfo<_dim, _type> \
           >::value                       \
         >::type* = nullptr>              \
void initField(SpheralTestField<T_DATA>& testField, size_t n)

#define TEST_RAW_DATA( _dim, _type)                     \
template<typename T_DATA,                               \
         typename std::enable_if<                       \
           std::is_same<                                \
             T_DATA,                                    \
             SpheralTypeInfo<_dim, _type>               \
           >::value                                     \
         >::type* = nullptr>                            \
void testSidreData(SpheralTestField<T_DATA>& testField, \
                   typename axom::sidre::Group *myFieldGroup)


#define INIT_ARITHMETIC_FIELD( _dim)                                             \
template<typename T_DATA,                                                        \
         typename T,                                                             \
         typename std::enable_if<std::is_arithmetic<T>::value>::type* = nullptr> \
void initField(Spheral::Field<typename T_DATA::DIM_T, T >& testField, size_t n)

#define TEST_ARITHMETIC_RAW_DATA( _dim)                                                                       \
template<typename T_DATA, typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type* = nullptr> \
void testSidreData(Spheral::Field<typename T_DATA::DIM_T, T >& testField,                                     \
                   typename axom::sidre::Group *myFieldGroup)


#define INIT_VECTOR_FIELD( _dim) \
template<typename T_DATA,        \
         typename T>             \
void initField(Spheral::Field<typename T_DATA::DIM_T, std::vector<T> >& testField, size_t n)

#define TEST_VECTOR_RAW_DATA( _dim)                                                    \
template<typename T_DATA, typename T>                                                  \
void testSidreData(Spheral::Field<typename T_DATA::DIM_T, std::vector<T> >& testField, \
                   typename axom::sidre::Group *myFieldGroup)


#define INIT_TUPLE_THREE_FIELD( _dim)  \
template<typename T_DATA,              \
         typename T>                   \
void initField(Spheral::Field<typename T_DATA::DIM_T, std::tuple<T, T, T> >& testField, size_t n)

#define TEST_TUPLE_THREE_RAW_DATA( _dim)                                                    \
template<typename T_DATA, typename T>                                                       \
void testSidreData(Spheral::Field<typename T_DATA::DIM_T, std::tuple<T, T, T> >& testField, \
                   typename axom::sidre::Group *myFieldGroup)


#define INIT_TUPLE_FOUR_FIELD( _dim)  \
template<typename T_DATA,             \
         typename T>                  \
void initField(Spheral::Field<typename T_DATA::DIM_T, std::tuple<T, T, T, T> >& testField, size_t n)

#define TEST_TUPLE_FOUR_RAW_DATA( _dim)                                                        \
template<typename T_DATA, typename T>                                                          \
void testSidreData(Spheral::Field<typename T_DATA::DIM_T, std::tuple<T, T, T, T> >& testField, \
                   typename axom::sidre::Group *myFieldGroup)


#define INIT_TUPLE_FIVE_FIELD( _dim)  \
template<typename T_DATA,             \
         typename T>                  \
void initField(Spheral::Field<typename T_DATA::DIM_T, std::tuple<T, T, T, T, T> >& testField, size_t n)

#define TEST_TUPLE_FIVE_RAW_DATA( _dim)                                                           \
template<typename T_DATA, typename T>                                                             \
void testSidreData(Spheral::Field<typename T_DATA::DIM_T, std::tuple<T, T, T, T, T> >& testField, \
                   typename axom::sidre::Group *myFieldGroup)

// ------------------------------------------------------------
// Dim 1 : T
// ------------------------------------------------------------
INIT_ARITHMETIC_FIELD(1)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = (T)i;
}

TEST_ARITHMETIC_RAW_DATA(1)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    T *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i], *rawSidreData);
  }
}

// ------------------------------------------------------------
// Dim 1 : String
// ------------------------------------------------------------
INIT_FIELD(1, std::string)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = "This is a test string: " + std::to_string(i) + "\n";
}

TEST_RAW_DATA(1, std::string)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    char *rawSidreData = myFieldGroup->getView(i)->getData();
    for (u_long j = 0; j < testField[i].size(); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}

// ------------------------------------------------------------
// Dim 1 : std::vector<T>
// ------------------------------------------------------------
INIT_VECTOR_FIELD(1)
{
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; j++)
      testField[i].emplace_back((T)i);
}

TEST_VECTOR_RAW_DATA(1)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    T *rawSidreData = myFieldGroup->getView(i)->getData();
    for (u_long j = 0; j < testField[i].size(); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}

// ------------------------------------------------------------
// Dim 1 : std::tuple<T, T, T>
// ------------------------------------------------------------
INIT_TUPLE_THREE_FIELD(1)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = std::make_tuple((T)i, (T)i, (T)i);
}

TEST_TUPLE_THREE_RAW_DATA(1)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    T *rawSidreData = myFieldGroup->getView(i)->getData();

    EXPECT_EQ(std::get<0>(testField[i]), rawSidreData[0]);
    EXPECT_EQ(std::get<1>(testField[i]), rawSidreData[1]);
    EXPECT_EQ(std::get<2>(testField[i]), rawSidreData[2]);
  }
}

// ------------------------------------------------------------
// Dim 1 : std::tuple<T, T, T, T>
// ------------------------------------------------------------
INIT_TUPLE_FOUR_FIELD(1)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = std::make_tuple((T)i, (T)i, (T)i, (T)i);
}

TEST_TUPLE_FOUR_RAW_DATA(1)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    T *rawSidreData = myFieldGroup->getView(i)->getData();

    EXPECT_EQ(std::get<0>(testField[i]), rawSidreData[0]);
    EXPECT_EQ(std::get<1>(testField[i]), rawSidreData[1]);
    EXPECT_EQ(std::get<2>(testField[i]), rawSidreData[2]);
    EXPECT_EQ(std::get<3>(testField[i]), rawSidreData[3]);
  }
}

// ------------------------------------------------------------
// Dim 1 : std::tuple<T, T, T, T, T>
// ------------------------------------------------------------
INIT_TUPLE_FIVE_FIELD(1)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = std::make_tuple((T)i, (T)i, (T)i, (T)i, (T)i);
}

TEST_TUPLE_FIVE_RAW_DATA(1)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    T *rawSidreData = myFieldGroup->getView(i)->getData();

    EXPECT_EQ(std::get<0>(testField[i]), rawSidreData[0]);
    EXPECT_EQ(std::get<1>(testField[i]), rawSidreData[1]);
    EXPECT_EQ(std::get<2>(testField[i]), rawSidreData[2]);
    EXPECT_EQ(std::get<3>(testField[i]), rawSidreData[3]);
    EXPECT_EQ(std::get<4>(testField[i]), rawSidreData[4]);
  }
}

// ------------------------------------------------------------
// Dim 1 : Vector
// ------------------------------------------------------------
INIT_FIELD(1, Spheral::Dim<1>::Vector)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<1>::Vector(i);
}

TEST_RAW_DATA(1, Spheral::Dim<1>::Vector)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].x(), rawSidreData[0]);
  }
}

// ------------------------------------------------------------
// Dim 1 : Vector3d
// ------------------------------------------------------------
INIT_FIELD(1, Spheral::Dim<1>::Vector3d)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<1>::Vector3d(i, i + 1, i + 2);
}

TEST_RAW_DATA(1, Spheral::Dim<1>::Vector3d)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].x(), rawSidreData[0]);
    EXPECT_EQ(testField[i].y(), rawSidreData[1]);
    EXPECT_EQ(testField[i].z(), rawSidreData[2]);
  }
}

// ------------------------------------------------------------
// Dim 1 : Tensor
// ------------------------------------------------------------
INIT_FIELD(1, Spheral::Dim<1>::Tensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<1>::Tensor(i);
}

TEST_RAW_DATA(1, Spheral::Dim<1>::Tensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].xx(), rawSidreData[0]);
  }
}

// ------------------------------------------------------------
// Dim 1 : SymTensor
// ------------------------------------------------------------
INIT_FIELD(1, Spheral::Dim<1>::SymTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<1>::SymTensor(i);
}

TEST_RAW_DATA(1, Spheral::Dim<1>::SymTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].xx(), rawSidreData[0]);
  }
}

// ------------------------------------------------------------
// Dim 1 : ThirdRankTensor
// ------------------------------------------------------------
INIT_FIELD(1, Spheral::Dim<1>::ThirdRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<1>::ThirdRankTensor(i);
}

TEST_RAW_DATA(1, Spheral::Dim<1>::ThirdRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i], rawSidreData[0]);
  }
}

// ------------------------------------------------------------
// Dim 1 : FourthRankTensor
// ------------------------------------------------------------
INIT_FIELD(1, Spheral::Dim<1>::FourthRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<1>::FourthRankTensor(i);
}

TEST_RAW_DATA(1, Spheral::Dim<1>::FourthRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i], rawSidreData[0]);
  }
}

// ------------------------------------------------------------
// Dim 1 : FifthRankTensor
// ------------------------------------------------------------
INIT_FIELD(1, Spheral::Dim<1>::FifthRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<1>::FifthRankTensor(i);
}

TEST_RAW_DATA(1, Spheral::Dim<1>::FifthRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i], rawSidreData[0]);
  }
}

// ------------------------------------------------------------
// Dim 2 : Vector
// ------------------------------------------------------------
INIT_FIELD(2, Spheral::Dim<2>::Vector)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<2>::Vector(i, i + 1);
}

TEST_RAW_DATA(2, Spheral::Dim<2>::Vector)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].x(), rawSidreData[0]);
    EXPECT_EQ(testField[i].y(), rawSidreData[1]);
  }
}

// ------------------------------------------------------------
// Dim 2 : Tensor
// ------------------------------------------------------------
INIT_FIELD(2, Spheral::Dim<2>::Tensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<2>::Tensor(i, i + 1, i + 2, i + 3);
}

TEST_RAW_DATA(2, Spheral::Dim<2>::Tensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].xx(), rawSidreData[0]);
    EXPECT_EQ(testField[i].xy(), rawSidreData[1]);
    EXPECT_EQ(testField[i].yx(), rawSidreData[2]);
    EXPECT_EQ(testField[i].yy(), rawSidreData[3]);
  }
}

// ------------------------------------------------------------
// Dim 2 : SymTensor
// ------------------------------------------------------------
INIT_FIELD(2, Spheral::Dim<2>::SymTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<2>::SymTensor(i, i + 1, i + 1, i);
}

TEST_RAW_DATA(2, Spheral::Dim<2>::SymTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].xx(), rawSidreData[0]);
    EXPECT_EQ(testField[i].xy(), rawSidreData[1]);
    EXPECT_EQ(testField[i].yy(), rawSidreData[2]);
  }
}

// ------------------------------------------------------------
// Dim 2 : ThirdRankTensor
// ------------------------------------------------------------
INIT_FIELD(2, Spheral::Dim<2>::ThirdRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<2>::ThirdRankTensor(i);
}

TEST_RAW_DATA(2, Spheral::Dim<2>::ThirdRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    for (int j = 0; j < Spheral::DataTypeTraits<Spheral::Dim<2>::ThirdRankTensor>::numElements(testField[0]); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}

// ------------------------------------------------------------
// Dim 2 : FourthRankTensor
// ------------------------------------------------------------
INIT_FIELD(2, Spheral::Dim<2>::FourthRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<2>::FourthRankTensor(i);
}

TEST_RAW_DATA(2, Spheral::Dim<2>::FourthRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    for (int j = 0; j < Spheral::DataTypeTraits<Spheral::Dim<2>::FourthRankTensor>::numElements(testField[0]); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}

// ------------------------------------------------------------
// Dim 2 : FifthRankTensor
// ------------------------------------------------------------
INIT_FIELD(2, Spheral::Dim<2>::FifthRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<2>::FifthRankTensor(i);
}

TEST_RAW_DATA(2, Spheral::Dim<2>::FifthRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    for (int j = 0; j < Spheral::DataTypeTraits<Spheral::Dim<2>::FifthRankTensor>::numElements(testField[0]); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}

// ------------------------------------------------------------
// Dim 3 : Vector
// ------------------------------------------------------------
INIT_FIELD(3, Spheral::Dim<3>::Vector)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<3>::Vector(i, i + 1, i + 2);
}

TEST_RAW_DATA(3, Spheral::Dim<3>::Vector)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].x(), rawSidreData[0]);
    EXPECT_EQ(testField[i].y(), rawSidreData[1]);
    EXPECT_EQ(testField[i].z(), rawSidreData[2]);
  }
}

// ------------------------------------------------------------
// Dim 3 : Tensor
// ------------------------------------------------------------
INIT_FIELD(3, Spheral::Dim<3>::Tensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<3>::Tensor(i, i + 1, i + 2,
                                           i, i + 1, i + 2,
                                           i, i + 1, i + 2);
}

TEST_RAW_DATA(3, Spheral::Dim<3>::Tensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].xx(), rawSidreData[0]);
    EXPECT_EQ(testField[i].xy(), rawSidreData[1]);
    EXPECT_EQ(testField[i].xz(), rawSidreData[2]);
    EXPECT_EQ(testField[i].yx(), rawSidreData[3]);
    EXPECT_EQ(testField[i].yy(), rawSidreData[4]);
    EXPECT_EQ(testField[i].yz(), rawSidreData[5]);
    EXPECT_EQ(testField[i].zx(), rawSidreData[6]);
    EXPECT_EQ(testField[i].zy(), rawSidreData[7]);
    EXPECT_EQ(testField[i].zz(), rawSidreData[8]);
  }
}

// ------------------------------------------------------------
// Dim 3 : SymTensor
// ------------------------------------------------------------
INIT_FIELD(3, Spheral::Dim<3>::SymTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<3>::SymTensor(i, i + 1, i + 2,
                                              i + 1, i, i + 1,
                                              i + 2, i + 1, i);
}

TEST_RAW_DATA(3, Spheral::Dim<3>::SymTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    EXPECT_EQ(testField[i].xx(), rawSidreData[0]);
    EXPECT_EQ(testField[i].xy(), rawSidreData[1]);
    EXPECT_EQ(testField[i].xz(), rawSidreData[2]);
    EXPECT_EQ(testField[i].yy(), rawSidreData[3]);
    EXPECT_EQ(testField[i].yz(), rawSidreData[4]);
    EXPECT_EQ(testField[i].zz(), rawSidreData[5]);
  }
}

// ------------------------------------------------------------
// Dim 3 : ThirdRankTensor
// ------------------------------------------------------------
INIT_FIELD(3, Spheral::Dim<3>::ThirdRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<3>::ThirdRankTensor(i);
}

TEST_RAW_DATA(3, Spheral::Dim<3>::ThirdRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    for (int j = 0; j < Spheral::DataTypeTraits<Spheral::Dim<3>::ThirdRankTensor>::numElements(testField[0]); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}

// ------------------------------------------------------------
// Dim 3 : FourthRankTensor
// ------------------------------------------------------------
INIT_FIELD(3, Spheral::Dim<3>::FourthRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<3>::FourthRankTensor(i);
}

TEST_RAW_DATA(3, Spheral::Dim<3>::FourthRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    for (int j = 0; j < Spheral::DataTypeTraits<Spheral::Dim<3>::FourthRankTensor>::numElements(testField[0]); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}

// ------------------------------------------------------------
// Dim 3 : FifthRankTensor
// ------------------------------------------------------------
INIT_FIELD(3, Spheral::Dim<3>::FifthRankTensor)
{
  for (size_t i = 0; i < n; ++i)
    testField[i] = Spheral::Dim<3>::FifthRankTensor(i);
}

TEST_RAW_DATA(3, Spheral::Dim<3>::FifthRankTensor)
{
  for (int i = 0; i < myFieldGroup->getNumViews(); ++i)
  {
    double *rawSidreData = myFieldGroup->getView(i)->getData();
    for (int j = 0; j < Spheral::DataTypeTraits<Spheral::Dim<3>::FifthRankTensor>::numElements(testField[0]); ++j)
      EXPECT_EQ(testField[i][j], rawSidreData[j]);
  }
}


using SpheralTypes = ::testing::Types<
  SpheralTypeInfo<1, char>,
  SpheralTypeInfo<1, int>,
  SpheralTypeInfo<1, size_t>,
  SpheralTypeInfo<1, uint32_t>,
  SpheralTypeInfo<1, uint64_t>,
  SpheralTypeInfo<1, float>,
  SpheralTypeInfo<1, double>,
  SpheralTypeInfo<1, std::string>,

  SpheralTypeInfo<1, std::vector<char>>,
  SpheralTypeInfo<1, std::vector<int>>,
  SpheralTypeInfo<1, std::vector<double>>,
  SpheralTypeInfo<1, std::vector<size_t>>,
  SpheralTypeInfo<1, std::vector<uint32_t>>,
  SpheralTypeInfo<1, std::vector<uint64_t>>,
  SpheralTypeInfo<1, std::vector<float>>,

  SpheralTypeInfo<1, std::tuple<char, char, char>>,
  SpheralTypeInfo<1, std::tuple<int, int, int>>,
  SpheralTypeInfo<1, std::tuple<size_t, size_t, size_t>>,
  SpheralTypeInfo<1, std::tuple<uint32_t, uint32_t, uint32_t>>,
  SpheralTypeInfo<1, std::tuple<uint64_t, uint64_t, uint64_t>>,
  SpheralTypeInfo<1, std::tuple<float, float, float>>,
  SpheralTypeInfo<1, std::tuple<double, double, double>>,

  SpheralTypeInfo<1, std::tuple<char, char, char, char>>,
  SpheralTypeInfo<1, std::tuple<int, int, int, int>>,
  SpheralTypeInfo<1, std::tuple<size_t, size_t, size_t, size_t>>,
  SpheralTypeInfo<1, std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>>,
  SpheralTypeInfo<1, std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>>,
  SpheralTypeInfo<1, std::tuple<float, float, float, float>>,
  SpheralTypeInfo<1, std::tuple<double, double, double, double>>,

  SpheralTypeInfo<1, std::tuple<char, char, char, char, char>>,
  SpheralTypeInfo<1, std::tuple<int, int, int, int, int>>,
  SpheralTypeInfo<1, std::tuple<size_t, size_t, size_t, size_t, size_t>>,
  SpheralTypeInfo<1, std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>>,
  SpheralTypeInfo<1, std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>>,
  SpheralTypeInfo<1, std::tuple<float, float, float, float, float>>,
  SpheralTypeInfo<1, std::tuple<double, double, double, double, double>>,
  
  SpheralTypeInfo<1, Spheral::Dim<1>::Vector>,
  SpheralTypeInfo<1, Spheral::Dim<1>::Vector3d>,
  SpheralTypeInfo<1, Spheral::Dim<1>::Tensor>,
  SpheralTypeInfo<1, Spheral::Dim<1>::SymTensor>,
  SpheralTypeInfo<1, Spheral::Dim<1>::ThirdRankTensor>,
  SpheralTypeInfo<1, Spheral::Dim<1>::FourthRankTensor>,
  SpheralTypeInfo<1, Spheral::Dim<1>::FifthRankTensor>,
  
  SpheralTypeInfo<2, Spheral::Dim<2>::Vector>,
  SpheralTypeInfo<2, Spheral::Dim<2>::Tensor>,
  SpheralTypeInfo<2, Spheral::Dim<2>::SymTensor>,
  SpheralTypeInfo<2, Spheral::Dim<2>::ThirdRankTensor>,
  SpheralTypeInfo<2, Spheral::Dim<2>::FourthRankTensor>,
  SpheralTypeInfo<2, Spheral::Dim<2>::FifthRankTensor>,
  
  SpheralTypeInfo<3, Spheral::Dim<3>::Vector>,
  SpheralTypeInfo<3, Spheral::Dim<3>::Tensor>,
  SpheralTypeInfo<3, Spheral::Dim<3>::SymTensor>,
  SpheralTypeInfo<3, Spheral::Dim<3>::ThirdRankTensor>,
  SpheralTypeInfo<3, Spheral::Dim<3>::FourthRankTensor>,
  SpheralTypeInfo<3, Spheral::Dim<3>::FifthRankTensor>
>;

template <typename T>
class SidreDataCollectionTestNew : public ::testing::Test {};
TYPED_TEST_SUITE(SidreDataCollectionTestNew, SpheralTypes);
