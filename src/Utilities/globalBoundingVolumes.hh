//---------------------------------Spheral++----------------------------------//
// globalBoundingVolumes
//
// Compute minimum bounding volumes (convex hull) for the nodes and their 
// sampling extents in the given DataBase.
// The two values computed are:
// 1.  nodeVolume -- the box containing the node positions.
// 2.  sampleVolume -- the box containing the node extents.
//
// Created by JMO, Sun Jan 31 19:53:36 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_globalBoundingVolumes__
#define __Spheral_globalBoundingVolumes__

#include <vector>
#include "Geometry/Dimension.hh"

// Forward declarations.
namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {

// Minimum bounding box for a Field of positions.
template<typename Dimension>
void
globalBoundingBox(const FieldSpace::Field<Dimension, typename Dimension::Vector>& positions,
                  typename Dimension::Vector& xmin,
                  typename Dimension::Vector& xmax,
                  const bool ghost = false);

// Minimum bounding box for a FieldList of positions.
template<typename Dimension>
void
globalBoundingBox(const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& positions,
                  typename Dimension::Vector& xmin,
                  typename Dimension::Vector& xmax,
                  const bool ghost = false);

// Minimum volume FacetedVolume for DataBases.
template<typename Dimension>
void 
globalBoundingVolumes(const DataBaseSpace::DataBase<Dimension>& dataBase,
                      typename Dimension::ConvexHull& nodeVolume,
                      typename Dimension::ConvexHull& sampleVolume);

}

#endif
