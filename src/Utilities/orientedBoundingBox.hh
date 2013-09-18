//---------------------------------Spheral++----------------------------------//
// orientedBoundingBox
//
// A collection of methods to compute oriented bounding boxes for collections
// of points, and test for containment in such boxes.
// There are two values computed:
// 1.  nodeBox -- the box containing the node positions.
// 2.  sampleBox -- the box containing the node extents.
//
// Created by JMO, Sun Jan 24 16:13:16 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_orientedBoundingBox__
#define __Spheral_orientedBoundingBox__

#include <vector>
#include "Geometry/Dimension.hh"
#include "spheralWildMagicConverters.hh"

// Forward declarations.
namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {

// Minimum volume boxes for DataBases.
template<typename Dimension>
void 
orientedBoundingBox(const DataBaseSpace::DataBase<Dimension>& dataBase,
                    typename Dimension::Box& nodeBox,
                    typename Dimension::Box& sampleBox);

// Minimum volume boxes for arrays of WMVectors.
Dim<1>::Box
orientedBoundingBox(const std::vector<Dim<1>::WMVector>& points);

Dim<2>::Box
orientedBoundingBox(const std::vector<Dim<2>::WMVector>& points);

Dim<3>::Box
orientedBoundingBox(const std::vector<Dim<3>::WMVector>& points);

// Minimum volume boxes for arrays of Vectors.
template<typename Dimension>
inline
typename Dimension::Box
orientedBoundingBox(const std::vector<typename Dimension::Vector>& points) {
  std::vector<typename Dimension::WMVector> WMpoints(points.size());
  std::transform(points.begin(), points.end(), WMpoints.begin(), convertVectorToWMVector<Dimension>);
  return orientedBoundingBox<Dimension>(WMpoints);
}

}

#endif
