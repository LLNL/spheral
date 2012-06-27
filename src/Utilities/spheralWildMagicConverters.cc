//---------------------------------Spheral++----------------------------------//
// spheralWildMagicConverters
//
// A collection of methods for converting assorted primitive types between 
// Spheral and WildMagic versions.
//
// Created by JMO, Mon Jan 25 16:28:37 PST 2010
//----------------------------------------------------------------------------//
#include <algorithm>
#include "spheralWildMagicConverters.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/DataBase.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Extract positions as WMVector from a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<typename Dimension::WMVector>
wildMagicPositions(const NodeList<Dimension>& nodes) {
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::WMVector WMVector;
  const Field<Dimension, Vector>& positions = nodes.positions();
  vector<WMVector> result;
  result.reserve(nodes.numNodes());
  transform(positions.begin(), positions.end(), back_inserter(result), convertVectorToWMVector<Dimension>);
  ENSURE(result.size() == nodes.numNodes());
  return result;
}

//------------------------------------------------------------------------------
// Extract positions as WMVectors from a DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<typename Dimension::WMVector>
wildMagicPositions(const DataBase<Dimension>& dataBase) {
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::WMVector WMVector;
  const int numNodes = dataBase.numNodes();
  const FieldList<Dimension, Vector> allpositions = dataBase.globalPosition();
  vector<WMVector> result;
  result.reserve(numNodes);
  for (typename FieldList<Dimension, Vector>::const_iterator itr = allpositions.begin();
       itr != allpositions.end();
       ++itr) {
    const Field<Dimension, Vector>& positions = **itr;
    transform(positions.begin(), positions.end(), back_inserter(result), convertVectorToWMVector<Dimension>);
  }
  ENSURE(result.size() == numNodes);
  return result;
}

//------------------------------------------------------------------------------
// Extract sampling boxes as WildMagic vectors from a NodeList (1-D).
//------------------------------------------------------------------------------
vector<Dim<1>::WMVector>
wildMagicSamplingPositions(const NodeList<Dim<1> >& nodes) {
  typedef Dim<1> Dimension;
  typedef Dimension::Vector Vector;
  typedef Dimension::WMVector WMVector;
  const Field<Dimension, Vector>& positions = nodes.positions();
  const Field<Dimension, Vector>& extents = nodes.neighbor().nodeExtentField();
  vector<WMVector> result;
  const int numNodes = nodes.numNodes();
  result.reserve(2*numNodes);
  for (int i = 0; i != numNodes; ++i) {
    const Vector& ri = positions(i);
    const Vector& extenti = extents(i);
    const double xi = ri.x();
    const double xext = extenti.x();
    result.push_back(WMVector(xi - xext));
    result.push_back(WMVector(xi + xext));
  }
  ENSURE(result.size() == 2*numNodes);
  return result;
}

//------------------------------------------------------------------------------
// Extract sampling boxes as WildMagic vectors from a NodeList (2-D).
//------------------------------------------------------------------------------
vector<Dim<2>::WMVector>
wildMagicSamplingPositions(const NodeList<Dim<2> >& nodes) {
  typedef Dim<2> Dimension;
  typedef Dimension::Vector Vector;
  typedef Dimension::WMVector WMVector;
  const Field<Dimension, Vector>& positions = nodes.positions();
  const Field<Dimension, Vector>& extents = nodes.neighbor().nodeExtentField();
  vector<WMVector> result;
  const int numNodes = nodes.numNodes();
  result.reserve(4*numNodes);
  for (int i = 0; i != numNodes; ++i) {
    const Vector& ri = positions(i);
    const Vector& extenti = extents(i);
    const double xi = ri.x();
    const double yi = ri.y();
    const double xext = extenti.x();
    const double yext = extenti.y();
    result.push_back(WMVector(xi - xext, yi - yext));
    result.push_back(WMVector(xi + xext, yi - yext));
    result.push_back(WMVector(xi - xext, yi + yext));
    result.push_back(WMVector(xi + xext, yi + yext));
  }
  ENSURE(result.size() == 4*numNodes);
  return result;
}

//------------------------------------------------------------------------------
// Extract sampling boxes as WildMagic vectors from a NodeList (3-D).
//------------------------------------------------------------------------------
vector<Dim<3>::WMVector >
wildMagicSamplingPositions(const NodeList<Dim<3> >& nodes) {
  typedef Dim<3> Dimension;
  typedef Dimension::Vector Vector;
  typedef Dimension::WMVector WMVector;
  const Field<Dimension, Vector>& positions = nodes.positions();
  const Field<Dimension, Vector>& extents = nodes.neighbor().nodeExtentField();
  vector<WMVector> result;
  const int numNodes = nodes.numNodes();
  result.reserve(8*numNodes);
  for (int i = 0; i != numNodes; ++i) {
    const Vector& ri = positions(i);
    const Vector& extenti = extents(i);
    const double xi = ri.x();
    const double yi = ri.y();
    const double zi = ri.z();
    const double xext = extenti.x();
    const double yext = extenti.y();
    const double zext = extenti.z();
    result.push_back(WMVector(xi - xext, yi - yext, zi - zext));
    result.push_back(WMVector(xi + xext, yi - yext, zi - zext));
    result.push_back(WMVector(xi - xext, yi + yext, zi - zext));
    result.push_back(WMVector(xi + xext, yi + yext, zi - zext));
    result.push_back(WMVector(xi - xext, yi - yext, zi + zext));
    result.push_back(WMVector(xi + xext, yi - yext, zi + zext));
    result.push_back(WMVector(xi - xext, yi + yext, zi + zext));
    result.push_back(WMVector(xi + xext, yi + yext, zi + zext));
  }
  ENSURE(result.size() == 8*numNodes);
  return result;
}

//------------------------------------------------------------------------------
// Extract sampling boxes as WMVectors from a DataBase (1-D).
//------------------------------------------------------------------------------
vector<Dim<1>::WMVector>
wildMagicSamplingPositions(const DataBase<Dim<1> >& dataBase) {
  typedef Dim<1> Dimension;
  typedef Dimension::Vector Vector;
  typedef Dimension::WMVector WMVector;
  const int numNodes = dataBase.numNodes();
  const FieldList<Dimension, Vector> allpositions = dataBase.globalPosition();
  const FieldList<Dimension, Vector> allextents = dataBase.globalNodeExtent();
  CHECK(allpositions.size() == allextents.size());
  vector<WMVector> result;
  result.reserve(2*numNodes);
  for (int ifield = 0; ifield != allpositions.size(); ++ifield) {
    const int n = allpositions[ifield]->numElements();
    for (int i = 0; i != n; ++i) {
      const Vector& ri = allpositions(ifield, i);
      const Vector& extenti = allextents(ifield, i);
      const double xi = ri.x();
      const double xext = extenti.x();
      result.push_back(WMVector(xi - xext));
      result.push_back(WMVector(xi + xext));
    }
  }
  ENSURE(result.size() == 2*numNodes);
  return result;
}

//------------------------------------------------------------------------------
// Extract sampling boxes as WMVectors from a DataBase (2-D).
//------------------------------------------------------------------------------
vector<Dim<2>::WMVector>
wildMagicSamplingPositions(const DataBase<Dim<2> >& dataBase) {
  typedef Dim<2> Dimension;
  typedef Dimension::Vector Vector;
  typedef Dimension::WMVector WMVector;
  const int numNodes = dataBase.numNodes();
  const FieldList<Dimension, Vector> allpositions = dataBase.globalPosition();
  const FieldList<Dimension, Vector> allextents = dataBase.globalNodeExtent();
  CHECK(allpositions.size() == allextents.size());
  vector<WMVector> result;
  result.reserve(4*numNodes);
  for (int ifield = 0; ifield != allpositions.size(); ++ifield) {
    const int n = allpositions[ifield]->numElements();
    for (int i = 0; i != n; ++i) {
      const Vector& ri = allpositions(ifield, i);
      const Vector& extenti = allextents(ifield, i);
      const double xi = ri.x();
      const double yi = ri.y();
      const double xext = extenti.x();
      const double yext = extenti.y();
      result.push_back(WMVector(xi - xext, yi - yext));
      result.push_back(WMVector(xi + xext, yi - yext));
      result.push_back(WMVector(xi - xext, yi + yext));
      result.push_back(WMVector(xi + xext, yi + yext));
    }
  }
  ENSURE(result.size() == 4*numNodes);
  return result;
}

//------------------------------------------------------------------------------
// Extract sampling boxes as WMVectors from a DataBase (3-D).
//------------------------------------------------------------------------------
vector<Dim<3>::WMVector>
wildMagicSamplingPositions(const DataBase<Dim<3> >& dataBase) {
  typedef Dim<3> Dimension;
  typedef Dimension::Vector Vector;
  typedef Dimension::WMVector WMVector;
  const int numNodes = dataBase.numNodes();
  const FieldList<Dimension, Vector> allpositions = dataBase.globalPosition();
  const FieldList<Dimension, Vector> allextents = dataBase.globalNodeExtent();
  CHECK(allpositions.size() == allextents.size());
  vector<WMVector> result;
  result.reserve(8*numNodes);
  for (int ifield = 0; ifield != allpositions.size(); ++ifield) {
    const int n = allpositions[ifield]->numElements();
    for (int i = 0; i != n; ++i) {
      const Vector& ri = allpositions(ifield, i);
      const Vector& extenti = allextents(ifield, i);
      const double xi = ri.x();
      const double yi = ri.y();
      const double zi = ri.z();
      const double xext = extenti.x();
      const double yext = extenti.y();
      const double zext = extenti.z();
      result.push_back(WMVector(xi - xext, yi - yext, zi - zext));
      result.push_back(WMVector(xi + xext, yi - yext, zi - zext));
      result.push_back(WMVector(xi - xext, yi + yext, zi - zext));
      result.push_back(WMVector(xi + xext, yi + yext, zi - zext));
      result.push_back(WMVector(xi - xext, yi - yext, zi + zext));
      result.push_back(WMVector(xi + xext, yi - yext, zi + zext));
      result.push_back(WMVector(xi - xext, yi + yext, zi + zext));
      result.push_back(WMVector(xi + xext, yi + yext, zi + zext));
    }
  }
  ENSURE(result.size() == 8*numNodes);
  return result;
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {

template vector<Dim<1>::WMVector> wildMagicPositions<Dim<1> >(const NodeList<Dim<1> >& nodes);
template vector<Dim<2>::WMVector> wildMagicPositions<Dim<2> >(const NodeList<Dim<2> >& nodes);
template vector<Dim<3>::WMVector> wildMagicPositions<Dim<3> >(const NodeList<Dim<3> >& nodes);

template vector<Dim<1>::WMVector> wildMagicPositions<Dim<1> >(const DataBase<Dim<1> >& dataBase);
template vector<Dim<2>::WMVector> wildMagicPositions<Dim<2> >(const DataBase<Dim<2> >& dataBase);
template vector<Dim<3>::WMVector> wildMagicPositions<Dim<3> >(const DataBase<Dim<3> >& dataBase);

}
