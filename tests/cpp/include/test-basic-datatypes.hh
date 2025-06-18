#ifndef SPHERAL_BASIC_DATATYPES_HH
#define SPHERAL_BASIC_DATATYPES_HH

#include "Utilities/DataTypeTraits.hh"
#include "config.hh"

// clang-format off
using TEST_FIELD_DATATYPES = camp::list<
  //bool
  char
  ,int
  ,size_t
  ,float
  ,double
  //,Spheral::Dim<3>::Vector
  //,Spheral::Dim<3>::Vector3d
  //,Spheral::Dim<1>::Tensor
  //Dim<1>::SymTensor,
  //Dim<1>::ThirdRankTensor,
  //Dim<1>::FourthRankTensor,
  //Dim<1>::FifthRankTensor
>;
// clang-format on

#endif // SPHERAL_BASIC_DATATYPES_HH
