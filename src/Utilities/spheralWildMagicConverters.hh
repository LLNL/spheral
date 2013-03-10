//---------------------------------Spheral++----------------------------------//
// spheralWildMagicConverters
//
// A collection of methods for converting assorted primitive types between 
// Spheral and WildMagic versions.
//
// Created by JMO, Mon Jan 25 16:28:37 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_spheralWildMagicConverters__
#define __Spheral_spheralWildMagicConverters__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

namespace NodeSpace {
  template<typename Dimension> class NodeList;
}
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}

//------------------------------------------------------------------------------
// Convert Spheral Vector -> WildMagic Vector.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::WMVector
convertVectorToWMVector(const typename Dimension::Vector& val);

// 1-D
template<>
inline
Dim<1>::WMVector
convertVectorToWMVector<Dim<1> >(const Dim<1>::Vector& val) {
  return val;
}

// 2-D
template<>
inline
Dim<2>::WMVector
convertVectorToWMVector<Dim<2> >(const Dim<2>::Vector& val) {
  return Dim<2>::WMVector(val.x(), val.y());
}

// 3-D
template<>
inline
Dim<3>::WMVector
convertVectorToWMVector<Dim<3> >(const Dim<3>::Vector& val) {
  return Dim<3>::WMVector(val.x(), val.y(), val.z());
}

//------------------------------------------------------------------------------
// Convert WildMagic Vector -> Spheral Vector.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
convertWMVectorToVector(const typename Dimension::WMVector& val);

template<>
inline
Dim<1>::Vector
convertWMVectorToVector<Dim<1> >(const Dim<1>::WMVector& val) {
  return val;
}

template<>
inline
Dim<2>::Vector
convertWMVectorToVector<Dim<2> >(const Dim<2>::WMVector& val) {
  return Dim<2>::Vector(val[0], val[1]);
}

template<>
inline
Dim<3>::Vector
convertWMVectorToVector<Dim<3> >(const Dim<3>::WMVector& val) {
  return Dim<3>::Vector(val[0], val[1], val[2]);
}

//------------------------------------------------------------------------------
// Extract all positions from a NodeList in WildMagic form.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename Dimension::WMVector>
wildMagicPositions(const NodeSpace::NodeList<Dimension>& nodes);

//------------------------------------------------------------------------------
// Extract all positions from a DataBase in WildMagic form.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename Dimension::WMVector>
wildMagicPositions(const DataBaseSpace::DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Extract all node sampling boxes as a set of positions in WildMagic form for
// a NodeList.
//------------------------------------------------------------------------------
std::vector<Dim<1>::WMVector>
wildMagicSamplingPositions(const NodeSpace::NodeList<Dim<1> >& nodes);

std::vector<Dim<2>::WMVector>
wildMagicSamplingPositions(const NodeSpace::NodeList<Dim<2> >& nodes);

std::vector<Dim<3>::WMVector>
wildMagicSamplingPositions(const NodeSpace::NodeList<Dim<3> >& nodes);

//------------------------------------------------------------------------------
// Extract all node sampling boxes as a set of positions in WildMagic form from
// a DataBase.
//------------------------------------------------------------------------------
std::vector<Dim<1>::WMVector>
wildMagicSamplingPositions(const DataBaseSpace::DataBase<Dim<1> >& dataBase);

std::vector<Dim<2>::WMVector>
wildMagicSamplingPositions(const DataBaseSpace::DataBase<Dim<2> >& dataBase);

std::vector<Dim<3>::WMVector>
wildMagicSamplingPositions(const DataBaseSpace::DataBase<Dim<3> >& dataBase);

}

#endif

