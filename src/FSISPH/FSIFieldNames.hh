//---------------------------------Spheral++----------------------------------//
// FSIFieldNames -- Fields names used in FSI package
//----------------------------------------------------------------------------//
#ifndef _Spheral_FSIFieldNames_
#define _Spheral_FSIFieldNames_

#include <string>

namespace Spheral {

struct FSIFieldNames {
  static const std::string pressureGradient;
  static const std::string specificThermalEnergyGradient;
  static const std::string interfaceFlags;
  static const std::string interfaceAreaVectors;
  static const std::string interfaceNormals;
  static const std::string interfaceAngles;
  static const std::string interfaceFraction;
  static const std::string interfaceSmoothness;
  static const std::string smoothedInterfaceNormals;
  static const std::string interfaceSmoothnessNormalization;
  //static const std::string inverseEquivalentDeviatoricStress;
};

}

#endif
