//---------------------------------Spheral++----------------------------------//
// DataTypeTraits.  A simple trait class to define the number of independent
// elements in the given DataType.
//
// Created by J. Michael Owen, Thu Oct 18 16:36:21 PDT 2001
//----------------------------------------------------------------------------//
#ifndef DataTypeTraits_HH
#define DataTypeTraits_HH

#include <stdint.h>
#include <vector>
#include "boost/tuple/tuple.hpp"
#include "Geometry/Dimension.hh"
#include "RegisterMPIDataTypes.hh"

#ifdef USE_MPI
extern "C" {
#include "mpi.h"
}
#endif

namespace Spheral {

template<typename DataType> struct DataTypeTraits;

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<char> {
  typedef char ElementType;
  static int numElements() { return 1; }
  static int zero() { return '\0'; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_CHAR; }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<int> {
  typedef int ElementType;
  static int numElements() { return 1; }
  static int zero() { return 0; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_INT; }
#endif
};

// //------------------------------------------------------------------------------
// template<>
// struct DataTypeTraits<unsigned> {
//   typedef unsigned ElementType;
//   static int numElements() { return 1; }
//   static unsigned zero() { return 0U; }
// #ifdef USE_MPI
//   static MPI_Datatype MpiDataType() { return MPI_UNSIGNED; }
// #endif
// };

// //------------------------------------------------------------------------------
// template<>
// struct DataTypeTraits<unsigned long> {
//   typedef unsigned long ElementType;
//   static int numElements() { return 1UL; }
//   static unsigned zero() { return 0UL; }
// #ifdef USE_MPI
//   static MPI_Datatype MpiDataType() { return MPI_UNSIGNED_LONG; }
// #endif
// };

// //------------------------------------------------------------------------------
// template<>
// struct DataTypeTraits<unsigned long long> {
//   typedef unsigned long long ElementType;
//   static int numElements() { return 1; }
//   static unsigned zero() { return 0ULL; }
// #ifdef USE_MPI
//   static MPI_Datatype MpiDataType() { return MPI_UNSIGNED_LONG_LONG; }
// #endif
// };

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<uint32_t> {
  typedef uint32_t ElementType;
  static int numElements() { return 1; }
  static uint32_t zero() { return 0UL; }
#ifdef USE_MPI
#ifdef MPI_UINT32_T
  static MPI_Datatype MpiDataType() { return MPI_UINT32_T; }
#else
  static MPI_Datatype MpiDataType() { return MPI_UNSIGNED; }
#endif
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<uint64_t> {
  typedef uint64_t ElementType;
  static int numElements() { return 1; }
  static uint64_t zero() { return 0ULL; }
#ifdef USE_MPI
#ifdef MPI_UINT64_T
  static MPI_Datatype MpiDataType() { return MPI_UINT64_T; }
#else
  static MPI_Datatype MpiDataType() { return MPI_UNSIGNED_LONG_LONG; }
#endif
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<float> {
  typedef float ElementType;
  static int numElements() { return 1; }
  static float zero() { return 0.0; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_FLOAT; }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<double> {
  typedef double ElementType;
  static int numElements() { return 1; }
  static double zero() { return 0.0; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_DOUBLE; }
#endif
};

//------------------------------------------------------------------------------
template<typename DataType>
struct DataTypeTraits<std::vector<DataType> > {
  typedef std::vector<DataType> ElementType;
  static int numElements() { return 0; }
  static std::vector<DataType> zero() { return std::vector<DataType>(); }
};

//------------------------------------------------------------------------------
template<typename DataType>
struct DataTypeTraits<boost::tuple<DataType, DataType, DataType> > {
  typedef boost::tuple<DataType, DataType, DataType> ElementType;
  static int numElements() { return 3; }
  static boost::tuple<DataType, DataType, DataType> zero() { return boost::make_tuple(DataType(), DataType(), DataType()); }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<DataType>::MpiDataType(); }
#endif
};

//------------------------------------------------------------------------------
template<typename DataType>
struct DataTypeTraits<boost::tuple<DataType, DataType, DataType, DataType> > {
  typedef boost::tuple<DataType, DataType, DataType, DataType> ElementType;
  static int numElements() { return 4; }
  static boost::tuple<DataType, DataType, DataType, DataType> zero() { return boost::make_tuple(DataType(), DataType(), DataType(), DataType()); }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<DataType>::MpiDataType(); }
#endif
};

//------------------------------------------------------------------------------
template<typename DataType>
struct DataTypeTraits<boost::tuple<DataType, DataType, DataType, DataType, DataType> > {
  typedef boost::tuple<DataType, DataType, DataType, DataType, DataType> ElementType;
  static int numElements() { return 5; }
  static boost::tuple<DataType, DataType, DataType, DataType, DataType> zero() { return boost::make_tuple(DataType(), DataType(), DataType(), DataType(), DataType()); }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<DataType>::MpiDataType(); }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<1>::Vector> {
  typedef double ElementType;
  static int numElements() { return 1; }
  static Dim<1>::Vector zero() { return Dim<1>::Vector::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::Vector3d> {
  typedef double ElementType;
  static int numElements() { return 3; }
  static Dim<1>::Vector3d zero() { return Dim<1>::Vector3d::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::Tensor> {
  typedef double ElementType;
  static int numElements() { return 1; }
  static Dim<1>::Tensor zero() { return Dim<1>::Tensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::SymTensor> {
  typedef double ElementType;
  static int numElements() { return 1; }
  static Dim<1>::SymTensor zero() { return Dim<1>::SymTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::ThirdRankTensor> {
  typedef double ElementType;
  static int numElements() { return Dim<1>::ThirdRankTensor::numElements; }
  static Dim<1>::ThirdRankTensor zero() { return Dim<1>::ThirdRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor1d; }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<2>::Vector> {
  typedef double ElementType;
  static int numElements() { return 2; }
  static Dim<2>::Vector zero() { return Dim<2>::Vector::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::Tensor> {
  typedef double ElementType;
  static int numElements() { return 4; }
  static Dim<2>::Tensor zero() { return Dim<2>::Tensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::SymTensor> {
  typedef double ElementType;
  static int numElements() { return 3; }
  static Dim<2>::SymTensor zero() { return Dim<2>::SymTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::ThirdRankTensor> {
  typedef double ElementType;
  static int numElements() { return Dim<2>::ThirdRankTensor::numElements; }
  static Dim<2>::ThirdRankTensor zero() { return Dim<2>::ThirdRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor2d; }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<3>::Vector> {
  typedef double ElementType;
  static int numElements() { return 3; }
  static Dim<3>::Vector zero() { return Dim<3>::Vector::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::Tensor> {
  typedef double ElementType;
  static int numElements() { return 9; }
  static Dim<3>::Tensor zero() { return Dim<3>::Tensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::SymTensor> {
  typedef double ElementType;
  static int numElements() { return 6; }
  static Dim<3>::SymTensor zero() { return Dim<3>::SymTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::ThirdRankTensor> {
  typedef double ElementType;
  static int numElements() { return Dim<3>::ThirdRankTensor::numElements; }
  static Dim<3>::ThirdRankTensor zero() { return Dim<3>::ThirdRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor3d; }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<1>::Box> {
  static Dim<1>::Box zero() { return Dim<1>::Box(Dim<1>::Vector(), 0.0); }
};

template<>
struct DataTypeTraits<Dim<2>::Box> {
  static Dim<2>::Box zero() { return Dim<2>::Box(); }
};

template<>
struct DataTypeTraits<Dim<3>::Box> {
  static Dim<3>::Box zero() { return Dim<3>::Box(); }
};

}

#else

namespace Spheral {
  template<typename DataType> struct DataTypeTraits;
}

#endif
