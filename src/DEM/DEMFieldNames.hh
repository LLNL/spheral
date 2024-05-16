//---------------------------------Spheral++----------------------------------//
// DEMFieldNames -- A collection of standard Field names for the DEM 
// physics package.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef _Spheral_DEMFieldNames_
#define _Spheral_DEMFieldNames_

#include <string>

namespace Spheral {

struct DEMFieldNames {
  static const std::string momentOfInertia;
  static const std::string particleRadius;
  static const std::string compositeParticleIndex;
  static const std::string angularVelocity;
  static const std::string uniqueIndices;
  static const std::string isActiveContact;
  static const std::string neighborIndices;
  static const std::string shearDisplacement;
  static const std::string rollingDisplacement;
  static const std::string torsionalDisplacement;
  static const std::string equilibriumOverlap;
  static const std::string maximumOverlap;
  static const std::string solidBoundaries;
  static const std::string solidBoundaryPolicy;
};

}

#endif
