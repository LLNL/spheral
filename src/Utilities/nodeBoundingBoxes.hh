//---------------------------------Spheral++----------------------------------//
// nodeBoundingBoxes
//
// Compute minimum volume bounding boxes for nodes and their extent.
//
// Created by JMO, Sun Jan 24 16:13:16 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_nodeBoundingBox__
#define __Spheral_nodeBoundingBox__

namespace Spheral {

// Forward declarations.
namespace NodeSpace {
  template<typename Dimension> class NodeList;
}
namespace FieldSpace {
  template<typename Dimension, typename Value> class Field;
  template<typename Dimension, typename Value> class FieldList;
}
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}

//------------------------------------------------------------------------------
// The bounding box for a position and H.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Box
boundingBox(const typename Dimension::Vector& xi,
            const typename Dimension::SymTensor& Hi,
            const typename Dimension::Scalar& kernelExtent);

//------------------------------------------------------------------------------
// The bounding boxes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldSpace::Field<Dimension, typename Dimension::Box>
nodeBoundingBoxes(const NodeSpace::NodeList<Dimension>& nodes);

//------------------------------------------------------------------------------
// The bounding boxes for all nodes in a DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldSpace::FieldList<Dimension, typename Dimension::Box>
nodeBoundingBoxes(const DataBaseSpace::DataBase<Dimension>& dataBase);

}

#include "nodeBoundingBoxesInline.hh"

#endif
