//------------------------------------------------------------------------------
// KernelIntegrationData
// 
// Integration information that is passed to integrals
//------------------------------------------------------------------------------
#ifndef __Spheral_KernelIntegrationData_hh__
#define __Spheral_KernelIntegrationData_hh__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension>
struct KernelIntegrationData {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  KernelIntegrationData() { }

  // Integration weight and ordinate in physical space
  // This could be either in the problem volume or on a surface
  double weight; // integration weight
  Vector ordinate; // integration point

  // Kernel value and dvalue for all points whose support includes integration point
  std::vector<Scalar> values; // values[k]
  std::vector<Vector> dvalues; // dvalues[k]
  
  // The index of the integration region, in case we want to use constant coefficients
  int index0;
  std::pair<int, int> nodeIndex0;

  // This converts between the integration index of the function and the flattened index
  // Could probably do without the nodeIndices if we changed the kernel evaluation
  std::vector<int> indices; // indices[integrationi] = locali
  std::vector<std::pair<int, int>> nodeIndices; // nodeIndices[integrationi] = (nodeListi, nodei)
  
  // The bilinear integrals are stored as integral[i][flatj]
  // Given the integration indices of i and j, return the flattened index for j
  std::vector<int> localIndex; // localIndex[integrationj + numLocalIndices * integrationj] = flatj

  // To interpolate with FieldLists, we need the volumes for each of the points
  std::vector<Scalar> volume; // volume[k]
  
  // If this is a surface integral, this is the normal for the integration surface
  Vector normal;

  // Surface integrals may depend on the surface being integrated
  // For linear integrals, the indexing is integral[i][flats]
  // For bilinear integrals, the indexing is integral[i][flats + numSurfacesi * flatj]
  int surfaceIndex0; // this is flats for the point that was used to create the integration cell
  std::vector<int> surfaceIndex; // surfaceIndex[integrationi] = flats

  // How many surfaces are there for the node at indices[i], used for bilinear/surface indexing
  std::vector<int> numSurfaces;

  // The current time, in case the coefficients need it
  Scalar time;
};

} // end namespace Spheral

#endif
