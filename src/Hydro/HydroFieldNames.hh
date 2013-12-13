//---------------------------------Spheral++----------------------------------//
// HydroFieldNames -- A collection of standard Field names for the hydro 
// physics package.
//
// Created by JMO, Sat Aug 28 21:16:14 2004
//----------------------------------------------------------------------------//
#ifndef _Spheral_HydroFieldNames_
#define _Spheral_HydroFieldNames_

#include <string>

namespace Spheral {

struct HydroFieldNames {
  static const std::string mass;
  static const std::string position;
  static const std::string velocity;
  static const std::string H;
  static const std::string work;
  static const std::string velocityGradient;
  static const std::string internalVelocityGradient;
  static const std::string massDensity;
  static const std::string specificThermalEnergy;
  static const std::string maxViscousPressure;
  static const std::string XSPHDeltaV;
  static const std::string XSPHWeightSum;
  static const std::string Hsmooth;
  static const std::string massFirstMoment;
  static const std::string massSecondMoment;
  static const std::string weightedNeighborSum;
  static const std::string pressure;
  static const std::string temperature;
  static const std::string soundSpeed;
  static const std::string pairAccelerations;
  static const std::string pairWork;
  static const std::string omegaGradh;
  static const std::string numberDensitySum;
  static const std::string timeStepMask;
  static const std::string A_CSPH;
  static const std::string B_CSPH;
  static const std::string C_CSPH;
  static const std::string D_CSPH;
  static const std::string gradA_CSPH;
  static const std::string gradB_CSPH;
  static const std::string volume;
  static const std::string linearMomentum;
  static const std::string totalEnergy;
  static const std::string mesh;
  static const std::string hourglassMask;
  static const std::string faceVelocity;
  static const std::string faceForce;
  static const std::string faceMass;
};

}

#else

namespace Spheral {
  string HydroFieldNames;
}

#endif
