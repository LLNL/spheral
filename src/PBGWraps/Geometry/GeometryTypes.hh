#ifndef __PBGWRAPS_GEOMETRYTYPES__
#define __PBGWRAPS_GEOMETRYTYPES__

#include "Python.h"
#include <vector>
#include <string>
#include <sstream>
#include "Geometry/Dimension.hh"
#include "Geometry/GeomVector.hh"
#include "Geometry/Geom3Vector.hh"
#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomThirdRankTensor.hh"
#include "Geometry/GeomFourthRankTensor.hh"
#include "Geometry/GeomFifthRankTensor.hh"
#include "Geometry/EigenStruct.hh"
#include "Geometry/computeEigenValues.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/GeomPolygon.hh"
#include "Geometry/GeomPolyhedron.hh"
#include "Geometry/invertRankNTensor.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerDoubleProduct.hh"
#include "Geometry/aggregateFacetedVolumes.hh"
#include "Utilities/DataTypeTraits.hh"

using namespace Spheral;

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef GeomVector<1> Vector1d;
typedef GeomVector<2> Vector2d;
typedef GeomVector<3> Vector3d;

typedef GeomTensor<1> Tensor1d;
typedef GeomTensor<2> Tensor2d;
typedef GeomTensor<3> Tensor3d;

typedef GeomSymmetricTensor<1> SymTensor1d;
typedef GeomSymmetricTensor<2> SymTensor2d;
typedef GeomSymmetricTensor<3> SymTensor3d;

typedef GeomThirdRankTensor<1> ThirdRankTensor1d;
typedef GeomThirdRankTensor<2> ThirdRankTensor2d;
typedef GeomThirdRankTensor<3> ThirdRankTensor3d;

typedef GeomFourthRankTensor<1> FourthRankTensor1d;
typedef GeomFourthRankTensor<2> FourthRankTensor2d;
typedef GeomFourthRankTensor<3> FourthRankTensor3d;

typedef GeomFifthRankTensor<1> FifthRankTensor1d;
typedef GeomFifthRankTensor<2> FifthRankTensor2d;
typedef GeomFifthRankTensor<3> FifthRankTensor3d;

typedef EigenStruct<1> EigenStruct1d;
typedef EigenStruct<2> EigenStruct2d;
typedef EigenStruct<3> EigenStruct3d;

typedef GeomPlane<Dim<1> > Plane1d;
typedef GeomPlane<Dim<2> > Plane2d;
typedef GeomPlane<Dim<3> > Plane3d;

typedef GeomFacet2d Facet2d;
typedef GeomPolygon Polygon;

typedef GeomFacet3d Facet3d;
typedef GeomPolyhedron Polyhedron;

}

typedef std::vector<Spheral::Facet2d> vector_of_Facet2d;
typedef std::vector<Spheral::Facet3d> vector_of_Facet3d;

typedef std::vector<Spheral::Box1d> vector_of_FacetedVolume1d;
typedef std::vector<Spheral::Polygon> vector_of_FacetedVolume2d;
typedef std::vector<Spheral::Polyhedron> vector_of_FacetedVolume3d;

typedef std::vector<std::vector<Spheral::Box1d> > vector_of_vector_of_FacetedVolume1d;
typedef std::vector<std::vector<Spheral::Polygon> > vector_of_vector_of_FacetedVolume2d;
typedef std::vector<std::vector<Spheral::Polyhedron> > vector_of_vector_of_FacetedVolume3d;

typedef std::vector<Spheral::Plane1d> vector_of_Plane1d;
typedef std::vector<Spheral::Plane2d> vector_of_Plane2d;
typedef std::vector<Spheral::Plane3d> vector_of_Plane3d;

//------------------------------------------------------------------------------
// Sequence methods for geometric types.
//------------------------------------------------------------------------------
namespace Spheral {

template<typename Type>
inline
unsigned
sizeGeomType(Type& self) {
  return DataTypeTraits<Type>::numElements(self);
}

template<typename Type>
inline
double
indexGeomType(Type& self,
              int index) {
  const int n = DataTypeTraits<Type>::numElements(self);
  if (index < 0) index += n;
  if (index < n) {
    return *(self.begin() + index);
  } else {
    PyErr_SetString(PyExc_IndexError, "Index out of range");
    return 0.0;
  }
}

template<typename Type>
inline
std::vector<double>
sliceGeomType(Type& self,
              int i1,
              int i2) {
  std::vector<double> result;
  const int n = DataTypeTraits<Type>::numElements(self);
  if (i1 < 0) i1 += n;
  if (i2 < 0) i2 += n;
  for (int i = i1; i < i2; ++i) result.push_back(indexGeomType(self, i));
  return result;
}

template<typename Type>
inline
int
assignToGeomType(Type& self, 
                 int index,
                 const double value) {
  const int n = DataTypeTraits<Type>::numElements(self);
  if (index < 0) index += n;
  if (index >= n) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return -1;
  } else {
    *(self.begin() + index) = value;
    return 0;
  }
}

template<typename Type>
inline
int
containsGeomType(Type& self,
                 const double value) {
  int result = 0;
  const int n = DataTypeTraits<Type>::numElements(self);
  int i = 0;
  while (result == 0 and i++ < n) {
    result = ((*(self.begin() + i) == value) ? 1 : 0);
  }
  return result;
}

//------------------------------------------------------------------------------
// Set the given element of a second rank tensor.
//------------------------------------------------------------------------------
template<typename TT>
inline
void
assignSecondRankTensorElement(TT& self,
                              const size_t i,
                              const size_t j,
                              const double val) {
  VERIFY(i < TT::nDimensions and j < TT::nDimensions);
  self(i,j) = val;
}

//------------------------------------------------------------------------------
// Set the given element of a third rank tensor.
//------------------------------------------------------------------------------
template<typename TT>
inline
void
assignThirdRankTensorElement(TT& self,
                             const size_t i,
                             const size_t j,
                             const size_t k,
                             const double val) {
  VERIFY(i < TT::nDimensions and j < TT::nDimensions and k < TT::nDimensions);
  self(i,j,k) = val;
}

//------------------------------------------------------------------------------
// Set the given element of a fourth rank tensor.
//------------------------------------------------------------------------------
template<typename TT>
inline
void
assignFourthRankTensorElement(TT& self,
                              const size_t i,
                              const size_t j,
                              const size_t k,
                              const size_t m,
                              const double val) {
  VERIFY(i < TT::nDimensions and j < TT::nDimensions and k < TT::nDimensions and m < TT::nDimensions);
  self(i,j,k,m) = val;
}

//------------------------------------------------------------------------------
// Set the given element of a fifth rank tensor.
//------------------------------------------------------------------------------
template<typename TT>
inline
void
assignFifthRankTensorElement(TT& self,
                              const size_t i,
                              const size_t j,
                              const size_t k,
                              const size_t m,
                              const size_t n,
                              const double val) {
  VERIFY(i < TT::nDimensions and j < TT::nDimensions and k < TT::nDimensions and m < TT::nDimensions and n < TT::nDimensions);
  self(i,j,k,m,n) = val;
}

//------------------------------------------------------------------------------
// Nice string representations (Vector)
//------------------------------------------------------------------------------
template<typename Vector>
inline
std::string
printReprVector(const Vector& val) {
  std::stringstream s;
  s << "Vector" << Vector::nDimensions << "d( ";
  for (size_t i = 0; i != Vector::nDimensions; ++i) {
    s << val(i) << " ";
  }
  s << ")";
  return s.str();
}

//------------------------------------------------------------------------------
// Nice string representations (Tensor)
//------------------------------------------------------------------------------
template<typename Tensor>
inline
std::string
printReprTensor(const Tensor& val) {
  std::stringstream s;
  s << "Tensor" << Tensor::nDimensions << "d(";
  for (size_t row = 0; row != Tensor::nDimensions; ++row) {
    s << "( ";
    for (size_t col = 0; col != Tensor::nDimensions; ++col) {
      s << val(row, col) << " ";
    }
    s << ")";
  }
  s << ")";
  return s.str();
}

//------------------------------------------------------------------------------
// Provide a function to construct a Geometry type (Vector, Tensor, etc.) from
// a sequence.
//------------------------------------------------------------------------------
template<typename GeomType>
GeomType*
constructGeomTypeFromSequence(PyObject* obj) {
  GeomType* result = new GeomType;
  if (PySequence_Check(obj) != 1) {
    PyErr_SetString(PyExc_ValueError, "Error: attempt to construct a GeomType from a non-sequence");
    return result;
  }
  const Py_ssize_t n = PySequence_Size(obj);
  if (n != DataTypeTraits<GeomType>::numElements(*result)) {
    PyErr_SetString(PyExc_ValueError, "Error: attempt to construct a GeomType from a sequence of incorrect size.");
    return result;
  }
  for (Py_ssize_t i = 0; i != n; ++i) {
    PyObject* numobj = PySequence_GetItem(obj, i);
    if (PyNumber_Check(numobj) != 1) {
      PyErr_SetString(PyExc_ValueError, "Error: attempt to construct a GeomType from a sequence containing non-numbers.");
      return result;
    }
    *(result->begin() + i) = PyFloat_AsDouble(numobj);
  }
  return result;
}

}

#endif
