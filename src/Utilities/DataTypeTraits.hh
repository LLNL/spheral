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
#include "Geometry/polyclipper.hh"
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
  static bool fixedSize() { return true; }
  static int numElements(const ElementType& x) { return 1; }
  static int zero() { return '\0'; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_CHAR; }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<int> {
  typedef int ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType& x) { return 1; }
  static int zero() { return 0; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_INT; }
#endif
};

// //------------------------------------------------------------------------------
// template<>
// struct DataTypeTraits<size_t> {
//   typedef size_t ElementType;
//   static bool fixedSize() { return true; }
//   static size_t numElements(const ElementType& x) { return 1; }
//   static size_t zero() { return 0U; }
// #ifdef USE_MPI
//   static MPI_Datatype MpiDataType() { return MPI_UNSIGNED; }
// #endif
// };

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<uint32_t> {
  typedef uint32_t ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType& x) { return 1; }
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
  static bool fixedSize() { return true; }
  static int numElements(const ElementType& x) { return 1; }
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
  static bool fixedSize() { return true; }
  static int numElements(const ElementType& x) { return 1; }
  static float zero() { return 0.0; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_FLOAT; }
#endif
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<double> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const ElementType& x) { return 1; }
  static double zero() { return 0.0; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return MPI_DOUBLE; }
#endif
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<std::vector<Value> > {
  typedef std::vector<Value> ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const std::vector<Value>& x) { return x.size(); }
  static std::vector<Value> zero() { return std::vector<Value>(); }
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<boost::tuple<Value, Value, Value> > {
  typedef Value ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const boost::tuple<Value, Value, Value>& x) { return 3; }
  static boost::tuple<Value, Value, Value> zero() { return boost::make_tuple(Value(), Value(), Value()); }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<Value>::MpiDataType(); }
#endif
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<boost::tuple<Value, Value, Value, Value> > {
  typedef boost::tuple<Value, Value, Value, Value> ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const boost::tuple<Value, Value, Value, Value>& x) { return 4; }
  static boost::tuple<Value, Value, Value, Value> zero() { return boost::make_tuple(Value(), Value(), Value(), Value()); }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<Value>::MpiDataType(); }
#endif
};

//------------------------------------------------------------------------------
template<typename Value>
struct DataTypeTraits<boost::tuple<Value, Value, Value, Value, Value> > {
  typedef boost::tuple<Value, Value, Value, Value, Value> ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const boost::tuple<Value, Value, Value, Value, Value>& x) { return 5; }
  static boost::tuple<Value, Value, Value, Value, Value> zero() { return boost::make_tuple(Value(), Value(), Value(), Value(), Value()); }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return DataTypeTraits<Value>::MpiDataType(); }
#endif
};

//------------------------------------------------------------------------------
template<typename Value1, typename Value2>
struct DataTypeTraits<std::pair<Value1, Value2> > {
  typedef std::pair<Value1, Value2> ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const std::pair<Value1, Value2>& x) { return 2; }
  static std::pair<Value1, Value2> zero() { return std::make_pair(Value1(), Value2()); }
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<1>::Vector> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::Vector& x) { return Dim<1>::Vector::numElements; }
  static Dim<1>::Vector zero() { return Dim<1>::Vector::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::Vector3d> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::Vector3d& x) { return Dim<1>::Vector3d::numElements; }
  static Dim<1>::Vector3d zero() { return Dim<1>::Vector3d::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::Tensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::Tensor& x) { return Dim<1>::Tensor::numElements; }
  static Dim<1>::Tensor zero() { return Dim<1>::Tensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::SymTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::SymTensor& x) { return Dim<1>::SymTensor::numElements; }
  static Dim<1>::SymTensor zero() { return Dim<1>::SymTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::ThirdRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::ThirdRankTensor& x) { return Dim<1>::ThirdRankTensor::numElements; }
  static Dim<1>::ThirdRankTensor zero() { return Dim<1>::ThirdRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::FourthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::FourthRankTensor& x) { return Dim<1>::FourthRankTensor::numElements; }
  static Dim<1>::FourthRankTensor zero() { return Dim<1>::FourthRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FourthRankTensor1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::FifthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<1>::FifthRankTensor& x) { return Dim<1>::FifthRankTensor::numElements; }
  static Dim<1>::FifthRankTensor zero() { return Dim<1>::FifthRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FifthRankTensor1d; }
#endif
};

template<>
struct DataTypeTraits<Dim<1>::FacetedVolume> {
  static Dim<1>::FacetedVolume zero() { return Dim<1>::FacetedVolume(); }
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<2>::Vector> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::Vector& x) { return Dim<2>::Vector::numElements; }
  static Dim<2>::Vector zero() { return Dim<2>::Vector::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::Tensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::Tensor& x) { return Dim<2>::Tensor::numElements; }
  static Dim<2>::Tensor zero() { return Dim<2>::Tensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::SymTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::SymTensor& x) { return Dim<2>::SymTensor::numElements; }
  static Dim<2>::SymTensor zero() { return Dim<2>::SymTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::ThirdRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::ThirdRankTensor& x) { return Dim<2>::ThirdRankTensor::numElements; }
  static Dim<2>::ThirdRankTensor zero() { return Dim<2>::ThirdRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::FourthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::FourthRankTensor& x) { return Dim<2>::FourthRankTensor::numElements; }
  static Dim<2>::FourthRankTensor zero() { return Dim<2>::FourthRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FourthRankTensor2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::FifthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<2>::FifthRankTensor& x) { return Dim<2>::FifthRankTensor::numElements; }
  static Dim<2>::FifthRankTensor zero() { return Dim<2>::FifthRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FifthRankTensor2d; }
#endif
};

template<>
struct DataTypeTraits<Dim<2>::FacetedVolume> {
  static Dim<2>::FacetedVolume zero() { return Dim<2>::FacetedVolume(); }
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<Dim<3>::Vector> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::Vector& x) { return Dim<3>::Vector::numElements; }
  static Dim<3>::Vector zero() { return Dim<3>::Vector::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Vector3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::Tensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::Tensor& x) { return Dim<3>::Tensor::numElements; }
  static Dim<3>::Tensor zero() { return Dim<3>::Tensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_Tensor3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::SymTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::SymTensor& x) { return Dim<3>::SymTensor::numElements; }
  static Dim<3>::SymTensor zero() { return Dim<3>::SymTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_SymTensor3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::ThirdRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::ThirdRankTensor& x) { return Dim<3>::ThirdRankTensor::numElements; }
  static Dim<3>::ThirdRankTensor zero() { return Dim<3>::ThirdRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_ThirdRankTensor3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::FourthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::FourthRankTensor& x) { return Dim<3>::FourthRankTensor::numElements; }
  static Dim<3>::FourthRankTensor zero() { return Dim<3>::FourthRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FourthRankTensor3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::FifthRankTensor> {
  typedef double ElementType;
  static bool fixedSize() { return true; }
  static int numElements(const Dim<3>::FifthRankTensor& x) { return Dim<3>::FifthRankTensor::numElements; }
  static Dim<3>::FifthRankTensor zero() { return Dim<3>::FifthRankTensor::zero; }
#ifdef USE_MPI
  static MPI_Datatype MpiDataType() { return RegisterMPIDataTypes::instance().MPI_FifthRankTensor3d; }
#endif
};

template<>
struct DataTypeTraits<Dim<3>::FacetedVolume> {
  static Dim<3>::FacetedVolume zero() { return Dim<3>::FacetedVolume(); }
};

//------------------------------------------------------------------------------
template<>
struct DataTypeTraits<PolyClipper::Vertex2d> {
  typedef PolyClipper::Vertex2d ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const ElementType& x) { return (DataTypeTraits<Dim<2>::Vector>::numElements(Dim<2>::Vector::zero) + 4); }
  static ElementType zero() { return PolyClipper::Vertex2d(); }
};

template<>
struct DataTypeTraits<PolyClipper::Vertex3d> {
  typedef PolyClipper::Vertex3d ElementType;
  static bool fixedSize() { return false; }
  static int numElements(const ElementType& x) { return (DataTypeTraits<Dim<3>::Vector>::numElements(Dim<3>::Vector::zero) +
                                                         x.neighbors.size() +
                                                         2); }
  static ElementType zero() { return PolyClipper::Vertex3d(); }
};

}

#else

namespace Spheral {
  template<typename DataType> struct DataTypeTraits;
}

#endif
