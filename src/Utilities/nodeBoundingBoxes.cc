//---------------------------------Spheral++----------------------------------//
// nodeBoundingBoxes
//
// Compute minimum volume bounding boxes for nodes and their extent.
//
// Created by JMO, Sun Jan 24 16:13:16 PST 2010
//----------------------------------------------------------------------------//
#include "nodeBoundingBoxes.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/DataBase.hh"

namespace Spheral {

using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// The bounding boxes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Box>
nodeBoundingBoxes(const NodeSpace::NodeList<Dimension>& nodes) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::Box Box;

  Field<Dimension, Box> result("NodeList bounding boxes", nodes);
  const Field<Dimension, Vector>& positions = nodes.positions();
  const Field<Dimension, SymTensor>& Hfield = nodes.Hfield();
  const Scalar kernelExtent = nodes.neighbor().kernelExtent();
  for (int i = 0; i != nodes.numNodes(); ++i) {
    result(i) = boundingBox<Dimension>(positions(i), Hfield(i), kernelExtent);
  }
  return result;
}


//------------------------------------------------------------------------------
// The bounding boxes for all nodes in a DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Box>
nodeBoundingBoxes(const DataBaseSpace::DataBase<Dimension>& dataBase) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::Box Box;

  FieldList<Dimension, Box> result = dataBase.newGlobalFieldList(Box(), "Bounding boxes");
  const FieldList<Dimension, Vector> positions = dataBase.globalPosition();
  const FieldList<Dimension, SymTensor> Hfield = dataBase.globalHfield();
  const Scalar kernelExtent = dataBase.maxKernelExtent();
  int nodeListi = 0;
  for (typename DataBase<Dimension>::ConstNodeListIterator itr = dataBase.nodeListBegin();
       itr != dataBase.nodeListEnd();
       ++itr, ++nodeListi) {
    const int numNodes = (**itr).numNodes();
    const Scalar kernelExtent = (**itr).neighbor().kernelExtent();
    for (int i = 0; i != numNodes; ++i) {
      result(nodeListi, i) = boundingBox<Dimension>(positions(nodeListi, i), Hfield(nodeListi, i), kernelExtent);
    }
  }
  return result;
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {

template Field<Dim<1>, Dim<1>::Box> nodeBoundingBoxes(const NodeList<Dim<1> >& nodes);
template Field<Dim<2>, Dim<2>::Box> nodeBoundingBoxes(const NodeList<Dim<2> >& nodes);
template Field<Dim<3>, Dim<3>::Box> nodeBoundingBoxes(const NodeList<Dim<3> >& nodes);

template FieldList<Dim<1>, Dim<1>::Box> nodeBoundingBoxes(const DataBase<Dim<1> >& dataBase);
template FieldList<Dim<2>, Dim<2>::Box> nodeBoundingBoxes(const DataBase<Dim<2> >& dataBase);
template FieldList<Dim<3>, Dim<3>::Box> nodeBoundingBoxes(const DataBase<Dim<3> >& dataBase);

}
