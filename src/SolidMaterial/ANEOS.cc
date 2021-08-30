//---------------------------------Spheral++----------------------------------//
// ANEOS -- An interface to the ANEOS equation of state from 
// Melosh amony others.  This is a C++ wrapper around calls to the underlying
// ANEOS fortran library.  The user must provide the ancillary file containing
// parameters for ANEOS in its expected format.  Also make sure to keep units
// consistent in this separate input file with what you give the C++ EOS here!
//
// Created by JMO, Tue Apr 23 14:55:28 PDT 2013
//----------------------------------------------------------------------------//
#include "ANEOS.hh"
#include "Field/Field.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <iostream>
using std::vector;
using std::string;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;


// Fortran baby!
extern "C" {
  void aneos_initialize_(char* filename, int* num, int* izetl);
  void call_aneos_(int* matnum, double* T, double* rho, 
                   double* P, double* E, double* S, double* CV, double* DPDT, double* DPDR, double* cs);
  double get_aneos_atomicweight_(int* matnum);
  double get_aneos_referencedensity_(int* matnum);
}

namespace Spheral {

// //------------------------------------------------------------------------------
// // Define an inline common method to handle calling ANEOS1.
// //------------------------------------------------------------------------------
// inline
// void call_ANEOS1(double& T, double& rho, 
//                  double& P, double& E, double& S, double& CV, double& DPDT, double& DPDR,
//                  int& mat) {
//   anesqt_.sqts[0] = std::sqrt(T);  // What the?  Some stuff is passed in in the common block
//   aneos1_(&T, &rho, 
//           &P, &E, &S, &CV, &DPDT, &DPDR,
//           &mat);
// }

//------------------------------------------------------------------------------
// A functor to compute eps(rho, T) for use building the interpolation table
//------------------------------------------------------------------------------
class epsFunc {
public:
  epsFunc(int matNum, double rhoConv, double Tconv, double Econv):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mEconv(Econv) {};
  double operator()(const Dim<2>::Vector& x) const {
    rho = x[0]/mRhoConv;
    T = x[1]/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return eps * mEconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mEconv;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ANEOS<Dimension>::
ANEOS(const int materialNumber,
      const unsigned numRhoVals,
      const unsigned numTvals,
      const double rhoMin,
      const double rhoMax,
      const double Tmin,
      const double Tmax,
      const PhysicalConstants& constants,
      const double externalPressure,
      const double minimumPressure,
      const double maximumPressure,
      const MaterialPressureMinType minPressureType):
  SolidEquationOfState<Dimension>(get_aneos_referencedensity_(const_cast<int*>(&materialNumber)),  // not in the right units yet!
                                  0.0,                                           // dummy etamin
                                  std::numeric_limits<double>::max(),            // dummy etamax
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minPressureType),
  mMaterialNumber(materialNumber),
  mNumRhoVals(numRhoVals),
  mNumTvals(numTvals),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax),
  mTmin(Tmin),
  mTmax(Tmax),
  mExternalPressure(externalPressure),
  mSTEvals(boost::extents[numRhoVals][numTvals]),
  mEpsInterp(),
  mANEOSunits(0.01,   // cm expressed as meters.
              0.001,  // g expressed in kg.
              1.0),   // sec in secs.
  mRhoConv(1.0),
  mTconv(1.0),
  mPconv(1.0),
  mEconv(1.0),
  mCVconv(1.0),
  mVelConv(1.0),
  mSconv(1.0),
  mAtomicWeight(0.0) {
  VERIFY2(numRhoVals > 1,
          "ANEOS ERROR : specify numRhoVals > 1");
  VERIFY2(numTvals > 1,
          "ANEOS ERROR : specify numTvals > 1");
  VERIFY2(rhoMin < rhoMax,
          "ANEOS ERROR : specify rhoMin < rhoMax");
  VERIFY2(Tmin < Tmax,
          "ANEOS ERROR : specify Tmin < Tmax");
  VERIFY2(Tmin > 0.0,
          "ANEOS ERROR : specify Tmin > 0.0");

  // // Convert temperature range to log space.
  // mTmin = log(mTmin);
  // mTmax = log(mTmax);
  
  // Look up the atomic weight.
  mAtomicWeight = get_aneos_atomicweight_(&mMaterialNumber);
  VERIFY2(mAtomicWeight > 0.0, 
          "ANEOS ERROR : bad atomic weight for material " << mMaterialNumber << " : " << mAtomicWeight);

  // Build our unit conversion factors.  After looking through the ANEOS source some it appears to me 
  // that they use mostly CGS units, except for temperatures which are in eV.
  const double lconv = mANEOSunits.unitLengthMeters() / constants.unitLengthMeters(),
               mconv = mANEOSunits.unitMassKg() / constants.unitMassKg(),
               tconv = mANEOSunits.unitTimeSec() / constants.unitTimeSec();
  mRhoConv = mconv/(lconv*lconv*lconv);
  mPconv = mconv/(lconv*tconv*tconv);
  mTconv = 1.160452e4; // eV/kB (CRC 2013)
  mEconv = FastMath::square(lconv/tconv);
  mCVconv = mEconv/mTconv;
  mVelConv = lconv/tconv;
  mSconv = mEconv/(mconv*mTconv);

  // Fix the reference density.
  this->referenceDensity(this->referenceDensity() * mRhoConv);

  // Build the biquadratic interpolation function for eps(rho, T)
  mEpsInterp.initialize(Dim<2>::Vector(mRhoMin/mRhoConv, mTmin/mRhoConv),
                        Dim<2>::Vector(mRhoMax/mRhoConv, mTmax/mRhoConv),
                        mNumRhoVals,
                        mNumTvals,
                        epsFunc(mMaterialNumber, mRhoConv, mTconv, mEconv));

  // Build our lookup table.
  const double drho = (mRhoMax - mRhoMin)/(mNumRhoVals - 1u);
  const double dT = (mTmax - mTmin)/(mNumTvals - 1u);
  CHECK(drho > 0.0);
  CHECK(dT > 0.0);
  double rhoi, Ti;
  for (unsigned i = 0u; i < mNumRhoVals; ++i) {
    rhoi = min(mRhoMin + i*drho, mRhoMax);
    for (unsigned j = 0; j < mNumTvals; ++j) {
      Ti = min(mTmin + j*dT, mTmax);
      mSTEvals[i][j] = mEpsInterp(Dim<2>::Vector(rhoi, Ti));
    }
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ANEOS<Dimension>::
~ANEOS() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != (int)Pressure.size(); ++i) {
    Pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != (int)temperature.size(); ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  for (int i = 0; i != (int)specificThermalEnergy.size(); ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  for (int i = 0; i != (int)specificHeat.size(); ++i) {
    specificHeat(i) = this->specificHeat(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != (int)soundSpeed.size(); ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  int n = massDensity.numElements();
  if (n > 0) {
    Field<Dimension, Scalar> T("temperature", gamma.nodeList()),
                            cv("cv", gamma.nodeList());
    this->setTemperature(T, massDensity, specificThermalEnergy);
    this->setSpecificHeat(cv, massDensity, specificThermalEnergy);
    Scalar nDen;
    for (int i = 0; i != n; ++i) {
      nDen = massDensity(i)/mAtomicWeight;
      gamma(i) = 1.0 + mConstants.molarGasConstant()*nDen*safeInvVar(cv(i));
    }
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != (int)bulkModulus.size(); ++i) {
    bulkModulus(i)=this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ANEOS<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  for (int i = 0; i != (int)entropy.size(); ++i) {
    entropy(i)=this->entropy(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  double Ti, rhoi, Pi, Ei, Si, CVi, DPDTi, DPDRi, csi;
  rhoi = max(mRhoMin, min(mRhoMax, massDensity)) / mRhoConv;
  Ti = this->temperature(massDensity, specificThermalEnergy) / mTconv;
  call_aneos_(const_cast<int*>(&mMaterialNumber), &Ti, &rhoi,
              &Pi, &Ei, &Si, &CVi, &DPDTi, &DPDRi, &csi);

  // That's it.
  Pi *= mPconv;
  return this->applyPressureLimits(Pi - mExternalPressure);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// This is an inverse lookup in the eps(rho, T) table we've precomputed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {

  // Figure out where we are in the table.
  // const double logeps = log(specificThermalEnergy);
  const double drho = (mRhoMax - mRhoMin)/(mNumRhoVals - 1u);
  const double dT = (mTmax - mTmin)/(mNumTvals - 1u);
  const unsigned irho0 = min(mNumRhoVals - 2u, unsigned(max(0.0, (massDensity - mRhoMin)/drho))),
                 irho2 = irho0 + 2u;
  const_slice_type rho0_slice = mSTEvals[boost::indices[irho0][range(0, mNumTvals)]];
  const unsigned iT0 = max(0, min(int(mNumTvals - 2u), 
                                  bisectSearch(rho0_slice.begin(), rho0_slice.end(), specificThermalEnergy))),
                 iT2 = iT0 + 2u;

  // If we're off the table, do the best we can...
  if (massDensity < mRhoMin) {
    if (specificThermalEnergy < mSTEvals[irho0][iT0]) return mSTEvals[irho0][iT0];
  } else if (massDensity > mRhoMax) {
    if (specificThermalEnergy > mSTEvals[irho2][iT2]) return mSTEvals[irho2][iT2];
  } else if (iT0 == 0u) {
    if (

  VERIFY2((specificThermalEnergy - mSTEvals[irho0][iT0])*(specificThermalEnergy - mSTEvals[irho0][iT2]) <= 0.0,
          "Bad bounding on eps: " << irho0 << " " << iT0 << " " << specificThermalEnergy << " " << mSTEvals[irho0][iT0] << " " << mSTEvals[irho0][iT2]);

  // There are two possible values for the temperature from the quadratic equation, so pick the one bounded
  // by our table position.
  const auto& coeffs = mEpsInterp.coeffs();
  const auto  i0 = mEpsInterp.lowerBound(massDensity + 0.5*drho,
                                          mTmin + (iT0 + 0.5)*dT);
  const auto  B0 = coeffs[i0] + coeffs[i0+1u]*massDensity + coeffs[i0+4u]*massDensity*massDensity - specificThermalEnergy;
  const auto  B1 = coeffs[i0+2u] + coeffs[i0+3u]*massDensity;
  const auto  x = std::sqrt(max(0.0, B1*B1 - 4.0*B0*coeffs[i0+5u]));
  auto  result = (-B1 + x)*safeInvVar(2.0*coeffs[i0+5u]);
  if (result < mTmin + iT0*dT or result > mTmin + iT2*dT) {
    result = -(B1 + x)*safeInvVar(2.0*coeffs[i0+5u]);
  }
  VERIFY(result >= mTmin + iT0*dT and result <= mTmin + iT2*dT);
  // cerr << "Looking up temperature : " << massDensity << " " << specificThermalEnergy << endl
  //      << "                         " << mRhoMin << " " << mRhoMax << " : " << mTmin << " " << mTmax << endl
  //      << "                         " << irho0 << " " << iT0 << endl
  //      << "                         " << mSTEvals[irho0][iT0] << " " << mSTEvals[irho0][iT1] << " " << mSTEvals[irho1][iT0] << " " << mSTEvals[irho1][iT1] << endl
  //      << "                         " << u << " " << t << endl
  //      << "                         " << num << " " << den << endl
  //      << "                         " << mTmin + (iT0 + t)*dT << endl;
  return result;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  double Ti, rhoi, Pi, Ei, Si, CVi, DPDTi, DPDRi, csi;
  rhoi = max(mRhoMin, min(mRhoMax, massDensity)) / mRhoConv;
  Ti = temperature / mTconv;
  call_aneos_(const_cast<int*>(&mMaterialNumber), &Ti, &rhoi,
              &Pi, &Ei, &Si, &CVi, &DPDTi, &DPDRi, &csi);
  return Ei * mEconv;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  double Ti, rhoi, Pi, Ei, Si, CVi, DPDTi, DPDRi, csi;
  rhoi = max(mRhoMin, min(mRhoMax, massDensity)) / mRhoConv;
  Ti = max(mTmin, min(mTmax, temperature)) / mTconv;
  call_aneos_(const_cast<int*>(&mMaterialNumber), &Ti, &rhoi,
              &Pi, &Ei, &Si, &CVi, &DPDTi, &DPDRi, &csi);
  const auto nDen = massDensity/mAtomicWeight;
  return max(1.0e-1*mConstants.molarGasConstant()*nDen, CVi * mCVconv);
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  double Ti, rhoi, Pi, Ei, Si, CVi, DPDTi, DPDRi, csi;
  rhoi = max(mRhoMin, min(mRhoMax, massDensity)) / mRhoConv;
  Ti = this->temperature(massDensity, specificThermalEnergy) / mTconv;
  call_aneos_(const_cast<int*>(&mMaterialNumber), &Ti, &rhoi,
              &Pi, &Ei, &Si, &CVi, &DPDTi, &DPDRi, &csi);
  return csi * mVelConv;
  // return sqrt(max(1e-10, DPDRi)) * mVelConv;
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  const double Ti = this->temperature(massDensity, specificThermalEnergy);
  const double cvi = this->specificHeat(massDensity, Ti);
  const double nDen = massDensity/mAtomicWeight;
  return 1.0 + mConstants.molarGasConstant()*nDen*safeInvVar(cvi);
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  double Ti, rhoi, Pi, Ei, Si, CVi, DPDTi, DPDRi, csi;
  rhoi = max(mRhoMin, min(mRhoMax, massDensity)) / mRhoConv;
  Ti = this->temperature(massDensity, specificThermalEnergy) / mTconv;
  call_aneos_(const_cast<int*>(&mMaterialNumber), &Ti, &rhoi,
              &Pi, &Ei, &Si, &CVi, &DPDTi, &DPDRi, &csi);
  return std::abs(rhoi * DPDRi * mPconv);
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  double Ti, rhoi, Pi, Ei, Si, CVi, DPDTi, DPDRi, csi;
  rhoi = max(mRhoMin, min(mRhoMax, massDensity)) / mRhoConv;
  Ti = this->temperature(massDensity, specificThermalEnergy) / mTconv;
  call_aneos_(const_cast<int*>(&mMaterialNumber), &Ti, &rhoi,
              &Pi, &Ei, &Si, &CVi, &DPDTi, &DPDRi, &csi);

  // That's it.
  Si *= mSconv;
  return Si;
}

//------------------------------------------------------------------------------
// valid
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ANEOS<Dimension>::
valid() const {
  return true;
}

//------------------------------------------------------------------------------
// Access the attributes.
//------------------------------------------------------------------------------
template<typename Dimension>
int
ANEOS<Dimension>::
materialNumber() const {
  return mMaterialNumber;
}

template<typename Dimension>
unsigned
ANEOS<Dimension>::
numRhoVals() const {
  return mNumRhoVals;
}

template<typename Dimension>
unsigned
ANEOS<Dimension>::
numTvals() const {
  return mNumTvals;
}

template<typename Dimension>
double
ANEOS<Dimension>::
rhoMin() const {
  return mRhoMin;
}

template<typename Dimension>
double
ANEOS<Dimension>::
rhoMax() const {
  return mRhoMax;
}

template<typename Dimension>
double
ANEOS<Dimension>::
Tmin() const {
  return mTmin;
}

template<typename Dimension>
double
ANEOS<Dimension>::
Tmax() const {
  return mTmax;
}

template<typename Dimension>
const boost::multi_array<double, 2>&
ANEOS<Dimension>::
specificThermalEnergyVals() const {
  return mSTEvals;
}

//------------------------------------------------------------------------------
// External pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
double
ANEOS<Dimension>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension>
void
ANEOS<Dimension>::
externalPressure(const double x) {
  mExternalPressure = x;
}

//------------------------------------------------------------------------------
// Atomic weight.
//------------------------------------------------------------------------------
template<typename Dimension>
double
ANEOS<Dimension>::
atomicWeight() const {
  return mAtomicWeight;
}

}

