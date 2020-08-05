//---------------------------------Spheral++----------------------------------//
// FacetedVolumeUtilities -- Compute the generic ancillary data for polygons
// and polyhedra.
// These methods are generic regardless of dimensionality, so we put them here
// for sharing.
//
// Created by JMO, Sun Mar 23 11:24:18 PDT 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_Geometry_computeAncillaryGeometry__
#define __Spheral_Geometry_computeAncillaryGeometry__

#include "Utilities/safeInv.hh"

#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <deque>

namespace Spheral {
namespace GeometryUtilities {

template<typename PolyType>
inline
void
computeAncillaryGeometry(const PolyType& poly,
                         std::vector<std::vector<unsigned> >& vertexFacetConnectivity,
                         std::vector<std::vector<unsigned> >& facetFacetConnectivity,
                         std::vector<typename PolyType::Vector>& vertexUnitNormals,
                         bool computeFacetFacetConnectivity) {

  using std::vector;

  typedef typename PolyType::Vector Vector;
  typedef typename PolyType::Facet Facet;

  const auto& vertices = poly.vertices();
  const auto& facets = poly.facets();
  const auto nverts = vertices.size();
  const auto nfacets = facets.size();

  // Find the facets surrounding each vertex.
  vertexFacetConnectivity = vector<vector<unsigned> >(nverts);
  for (unsigned fi = 0; fi != nfacets; ++fi) {
    const Facet& facet = facets[fi];
    const vector<unsigned>& ipoints = facet.ipoints();
    for (const unsigned vi: ipoints) {
      vertexFacetConnectivity[vi].push_back(fi);
    }
  }
  for (auto i = 0u; i < nverts; ++i) {
    std::sort(vertexFacetConnectivity[i].begin(), vertexFacetConnectivity[i].end());
    vertexFacetConnectivity[i].erase(std::unique(vertexFacetConnectivity[i].begin(), vertexFacetConnectivity[i].end()), vertexFacetConnectivity[i].end());
  }
  CHECK(vertexFacetConnectivity.size() == nverts);

  // Construct the facet->facet connectivity.
  if (computeFacetFacetConnectivity) {
    facetFacetConnectivity = vector<vector<unsigned> >(nfacets);
    for (const auto& vertexFacets: vertexFacetConnectivity) {
      for (const auto fi: vertexFacets) {
        std::copy(vertexFacets.begin(), vertexFacets.end(), std::back_inserter(facetFacetConnectivity[fi]));
      }
    }
    for (auto i = 0u; i < nfacets; ++i) {
      std::sort(facetFacetConnectivity[i].begin(), facetFacetConnectivity[i].end());
      facetFacetConnectivity[i].erase(std::unique(facetFacetConnectivity[i].begin(), facetFacetConnectivity[i].end()), facetFacetConnectivity[i].end());
      CHECK(find(facetFacetConnectivity[i].begin(), facetFacetConnectivity[i].end(), i) != facetFacetConnectivity[i].end());
    }
    CHECK(facetFacetConnectivity.size() == nfacets);
  } else {
    facetFacetConnectivity.clear();
  }

  // Find the normals to the surface at each vertex.
  vertexUnitNormals = vector<Vector>(nverts);
  for (auto i = 0u; i < nverts; ++i) {
    for (const auto fi: vertexFacetConnectivity[i]) {
      vertexUnitNormals[i] += facets[fi].normal()*safeInvVar(facets[fi].area());
    }
    vertexUnitNormals[i] = vertexUnitNormals[i].unitVector();
  }
  CHECK(vertexUnitNormals.size() == nverts);
}

}
}

#endif
