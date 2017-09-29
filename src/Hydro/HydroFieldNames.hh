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
  static const std::string hydroAcceleration;
  static const std::string massDensity;
  static const std::string normalization;
  static const std::string specificThermalEnergy;
  static const std::string maxViscousPressure;
  static const std::string effectiveViscousPressure;
  static const std::string massDensityCorrection;
  static const std::string viscousWork;
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
  static const std::string gamma;
  static const std::string entropy;
  static const std::string PSPHcorrection;
  static const std::string omegaGradh;
  static const std::string numberDensitySum;
  static const std::string timeStepMask;
  static const std::string m0_CRKSPH;
  static const std::string m1_CRKSPH;
  static const std::string m2_CRKSPH;
  static const std::string m3_CRKSPH;
  static const std::string m4_CRKSPH;
  static const std::string gradM0_CRKSPH;
  static const std::string gradM1_CRKSPH;
  static const std::string gradM2_CRKSPH;
  static const std::string gradM3_CRKSPH;
  static const std::string gradM4_CRKSPH;
  static const std::string A0_CRKSPH;
  static const std::string A_CRKSPH;
  static const std::string B_CRKSPH;
  static const std::string C_CRKSPH;
  static const std::string gradA0_CRKSPH;
  static const std::string gradA_CRKSPH;
  static const std::string gradB_CRKSPH;
  static const std::string gradC_CRKSPH;
  static const std::string surfacePoint;
  static const std::string voidPoint;
  static const std::string etaVoidPoints;
  static const std::string M_SPHCorrection;
  static const std::string volume;
  static const std::string linearMomentum;
  static const std::string totalEnergy;
  static const std::string mesh;
  static const std::string hourglassMask;
  static const std::string faceVelocity;
  static const std::string faceForce;
  static const std::string faceMass;
  static const std::string polyvols;
  static const std::string massDensityGradient;
  static const std::string ArtificialViscousClMultiplier;
  static const std::string ArtificialViscousCqMultiplier;
};

}

#else

namespace Spheral {
  string HydroFieldNames;
}

#endif
