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
#include "Utilities/bisectRoot.hh"
#include "Utilities/DBC.hh"

#include "boost/multi_array.hpp"
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

namespace { // anonymous

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
// A functor to compute T(rho, eps) for use building the interpolation table
//------------------------------------------------------------------------------
class Tfunc {
public:
  Tfunc(double Tmin, 
        double Tmax,
        const BiQuadraticInterpolator& epsInterp):
    mTmin(Tmin),
    mTmax(Tmax),
    mEpsInterp(epsInterp) {}

  double operator()(const Dim<2>::Vector& x) const {
    // Check if we're in bounds.
    const auto rho = x[0], eps = x[1];
    const auto epsMin = mEpsInterp(rho, mTmin);
    const auto epsMax = mEpsInterp(rho, mTmax);
    if (eps < epsMin) {
      const auto cvInv = mTmin/epsMin;
      return max(0.0, mTmin + (eps - epsMin)*cvInv);
    } else if (eps > epsMax) {
      const auto cvInv = mTmax/epsMax;
      return max(0.0, mTmax + (eps - epsMax)*cvInv);
    } else {
      return bisectRoot(Trho_func(rho, eps, mEpsInterp),
                        mTmin, mTmax,
                        1.0e-15, 1.0e-15, 200u);
    }
  }

private:
  double mTmin, mTmax;
  const BiQuadraticInterpolator& mEpsInterp;

  // We need to make a single argument functor for eps(T) given a fixed rho
  class Trho_func {
    double mrho, meps;
    const BiQuadraticInterpolator& mEpsInterp;
  public:
    Trho_func(const double rho,
              const double eps,
              const BiQuadraticInterpolator& epsInterp):
      mrho(rho),
      meps(eps),
      mEpsInterp(epsInterp) {}
    double operator()(const double T) const { return mEpsInterp(mrho, T) - meps; }
  };
};

// //------------------------------------------------------------------------------
// // T(rho, eps) with extrapolation
// // Note this is only extrapolation in T, not rho
// //------------------------------------------------------------------------------
// class TfuncExtra {
// public:
//   TfuncExtra(const BiQuadraticInterpolator& Tinterp,
//              const BiQuadraticInterpolator& cvInterp,):
//     mTinterp(Tinterp),
//     mcvInterp(cvInterp) {}
//   double operator()(const Dim<2>::Vector& x) const {
//     rho = x[0];
//     eps = x[1];
//     if (eps < 

//      = mTinterp(x[0], x[1])/mTconv;
//     call_aneos_(&mMatNum, &T, &rho,
//                 &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
//     return P * mPconv;
//   }
// private:
//   mutable int mMatNum;
//   double mRhoConv, mTconv, mPconv;
//   const BiQuadraticInterpolator& mTinterp;
//   mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
// };

//------------------------------------------------------------------------------
// P(rho, eps)
//------------------------------------------------------------------------------
class Pfunc {
public:
  Pfunc(int matNum, double rhoConv, double Tconv, double Pconv,
        const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mTinterp(Tinterp) {}
  double operator()(const Dim<2>::Vector& x) const {
    rho = x[0]/mRhoConv;
    T = mTinterp(x[0], x[1])/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return P * mPconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const BiQuadraticInterpolator& mTinterp;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// cV(rho, eps)
//------------------------------------------------------------------------------
class cVfunc {
public:
  cVfunc(int matNum, double rhoConv, double Tconv, double cVconv):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mCVconv(cVconv) {}
  double operator()(const Dim<2>::Vector& x) const {
    rho = x[0]/mRhoConv;
    T = x[1]/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return cV * mCVconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mCVconv;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// cs(rho, eps)
//------------------------------------------------------------------------------
class csfunc {
public:
  csfunc(int matNum, double rhoConv, double Tconv, double velConv,
         const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mVelConv(velConv),
    mTinterp(Tinterp) {}
  double operator()(const Dim<2>::Vector& x) const {
    rho = x[0]/mRhoConv;
    T = mTinterp(x[0], x[1])/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return cs * mVelConv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mVelConv;
  const BiQuadraticInterpolator& mTinterp;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// DPDrho(rho, eps)
//------------------------------------------------------------------------------
class Kfunc {
public:
  Kfunc(int matNum, double rhoConv, double Tconv, double Pconv,
        const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mTinterp(Tinterp) {}
  double operator()(const Dim<2>::Vector& x) const {
    rho = x[0]/mRhoConv;
    T = mTinterp(x[0], x[1])/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return std::abs(rho * DPDR * mPconv);
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const BiQuadraticInterpolator& mTinterp;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// s(rho, eps) (entropy)
//------------------------------------------------------------------------------
class sfunc {
public:
  sfunc(int matNum, double rhoConv, double Tconv, double Sconv,
        const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mSconv(Sconv),
    mTinterp(Tinterp) {}
  double operator()(const Dim<2>::Vector& x) const {
    rho = x[0]/mRhoConv;
    T = mTinterp(x[0], x[1])/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return S * mSconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mSconv;
  const BiQuadraticInterpolator& mTinterp;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

} // anonymous

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
  mEpsMin(std::numeric_limits<double>::max()),
  mEpsMax(std::numeric_limits<double>::min()),
  mExternalPressure(externalPressure),
  mEpsInterp(),
  mTinterp(),
  mPinterp(),
  mCVinterp(),
  mCSinterp(),
  mKinterp(),
  mSinterp(),
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

  // Scan the boundaries of the (rho,T) space for the min/max specific thermal energy (eps)
  const auto Feps = epsFunc(mMaterialNumber, mRhoConv, mTconv, mEconv);
  const auto drho = (mRhoMax - mRhoMin)/(mNumRhoVals - 1u);
  const auto dT = (mTmax - mTmin)/(mNumTvals - 1u);
  CHECK(drho > 0.0);
  CHECK(dT > 0.0);
  // [rhoMin:rhoMax], Tmin
  for (auto i = 0u; i < mNumRhoVals; ++i) {
    const auto epsi = Feps(Dim<2>::Vector(mRhoMin + i*drho, mTmin));
    mEpsMin = min(mEpsMin, epsi);
    mEpsMax = max(mEpsMax, epsi);
  }
  // [rhoMin:rhoMax], Tmax
  for (auto i = 0u; i < mNumRhoVals; ++i) {
    const auto epsi = Feps(Dim<2>::Vector(mRhoMin + i*drho, mTmax));
    mEpsMin = min(mEpsMin, epsi);
    mEpsMax = max(mEpsMax, epsi);
  }
  // rhoMin, [Tmin:Tmax]
  for (auto i = 0u; i < mNumTvals; ++i) {
    const auto epsi = Feps(Dim<2>::Vector(mRhoMin, mTmin + i*dT));
    mEpsMin = min(mEpsMin, epsi);
    mEpsMax = max(mEpsMax, epsi);
  }
  // rhoMax, [Tmin:Tmax]
  for (auto i = 0u; i < mNumTvals; ++i) {
    const auto epsi = Feps(Dim<2>::Vector(mRhoMax, mTmin + i*dT));
    mEpsMin = min(mEpsMin, epsi);
    mEpsMax = max(mEpsMax, epsi);
  }

  // Build the biquadratic interpolation function for eps(rho, T)
  mEpsInterp.initialize(Dim<2>::Vector(mRhoMin, mTmin),
                        Dim<2>::Vector(mRhoMax, mTmax),
                        mNumRhoVals, mNumTvals, Feps);

  // Now the hard inversion method for looking up T(rho, eps)
  const auto Ftemp = Tfunc(mTmin, mTmax, mEpsInterp);
  mTinterp.initialize(Dim<2>::Vector(mRhoMin, mEpsMin),
                      Dim<2>::Vector(mRhoMax, mEpsMax),
                      mNumRhoVals, mNumTvals, Ftemp);

  // And finally the interpolators for most of our derived quantities
  const auto Fpres = Pfunc(mMaterialNumber, mRhoConv, mTconv, mPconv, mTinterp);
  mPinterp.initialize(Dim<2>::Vector(mRhoMin, mEpsMin),
                      Dim<2>::Vector(mRhoMax, mEpsMax),
                      mNumRhoVals, mNumTvals, Fpres);

  const auto Fcv = cVfunc(mMaterialNumber, mRhoConv, mTconv, mCVconv);
  mCVinterp.initialize(Dim<2>::Vector(mRhoMin, mTmin),
                       Dim<2>::Vector(mRhoMax, mTmax),
                       mNumRhoVals, mNumTvals, Fcv);

  const auto Fcs = csfunc(mMaterialNumber, mRhoConv, mTconv, mVelConv, mTinterp);
  mCSinterp.initialize(Dim<2>::Vector(mRhoMin, mEpsMin),
                       Dim<2>::Vector(mRhoMax, mEpsMax),
                       mNumRhoVals, mNumTvals, Fcs);

  const auto FK = Kfunc(mMaterialNumber, mRhoConv, mTconv, mPconv, mTinterp);
  mKinterp.initialize(Dim<2>::Vector(mRhoMin, mEpsMin),
                      Dim<2>::Vector(mRhoMax, mEpsMax),
                      mNumRhoVals, mNumTvals, FK);


  const auto Fs = sfunc(mMaterialNumber, mRhoConv, mTconv, mSconv, mTinterp);
  mSinterp.initialize(Dim<2>::Vector(mRhoMin, mEpsMin),
                      Dim<2>::Vector(mRhoMax, mEpsMax),
                      mNumRhoVals, mNumTvals, Fs);

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
  return this->applyPressureLimits(mPinterp(massDensity, specificThermalEnergy) - mExternalPressure);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  return mTinterp(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  return mEpsInterp(massDensity, temperature);
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  return mCVinterp(massDensity, temperature);
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  return mCSinterp(massDensity, specificThermalEnergy);
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
  return mKinterp(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  return mSinterp(massDensity, specificThermalEnergy);
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

