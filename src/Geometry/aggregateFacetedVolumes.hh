//---------------------------------Spheral++----------------------------------//
// aggregateFacetedVolumes
//
// Weld together FacetedVolumes into a single object.  Not really a union since
// we make no effort to remove redundancies or overlap.
// Only works with Polygons and Polyhedra.
//
// Created by JMO, Wed Mar  1 14:30:39 PST 2017
//----------------------------------------------------------------------------//
#ifndef __Spheral_aggregateFacetedVolumes__
#define __Spheral_aggregateFacetedVolumes__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension>
inline
typename Dimension::FacetedVolume 
aggregateFacetedVolumes(const std::vector<typename Dimension::FacetedVolume>& loops) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  std::vector<Vector> vertices;
  std::vector<std::vector<unsigned>> indices;
  for (const auto& loop: loops) {
    const auto& loopVerts = loop.vertices();
    const auto loopFacetVertices = loop.facetVertices();
    std::copy(loopVerts.begin(), loopVerts.end(), std::back_inserter(vertices));
    std::copy(loopFacetVertices.begin(), loopFacetVertices.end(), std::back_inserter(indices));
  }
  return FacetedVolume(vertices, indices);
}

}

#endif
