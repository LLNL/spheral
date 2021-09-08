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
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/bisectRoot.hh"
#include "Utilities/DBC.hh"

#include "boost/multi_array.hpp"
#include <iostream>
#include <ctime>
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
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = y/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return eps * mEconv;
  }

public:
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
        const QuadraticInterpolator& epsMinInterp,
        const BiQuadraticInterpolator& epsInterp,
        const epsFunc& Feps,
        const bool verbose = false):
    mTmin(Tmin),
    mTmax(Tmax),
    mEpsMinInterp(epsMinInterp),
    mEpsInterp(epsInterp),
    mFeps(Feps),
    mVerbose(verbose) {}

  double operator()(const double rho, const double eps) const {
    // Find the eps offset for this density
    const auto eps0 = mEpsMinInterp(rho);
    const auto FT = Trho_func(rho, eps + eps0, mEpsInterp);
    // cerr << " **> " << rho << " " << eps << " " << eps0 << " " << (eps + eps0) << " " << FT(mTmin) << " " << FT(mTmax) << endl;
    return bisectRoot(Trho_func(rho, eps + eps0, mEpsInterp),
                      mTmin, mTmax,
                      1.0e-15, 200u);
  }

private:
  double mTmin, mTmax;
  const QuadraticInterpolator& mEpsMinInterp;
  const BiQuadraticInterpolator& mEpsInterp;
  const epsFunc& mFeps;
  bool mVerbose;

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

//------------------------------------------------------------------------------
// P(rho, eps)
//------------------------------------------------------------------------------
class Pfunc {
public:
  Pfunc(int matNum, double rhoConv, double Tconv, double Pconv,
        const QuadraticInterpolator& epsMinInterp,
        const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mEpsMinInterp(epsMinInterp),
    mTinterp(Tinterp) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    const auto eps0 = mEpsMinInterp(x);
    // cerr << " *** mTinterp(" << x << ") = " << mTinterp(x[0], x[1]) << endl;
    T = mTinterp(x, y - eps0)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return P * mPconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const QuadraticInterpolator& mEpsMinInterp;
  const BiQuadraticInterpolator& mTinterp;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// cV(rho, T)
//------------------------------------------------------------------------------
class cVfunc {
public:
  cVfunc(int matNum, double rhoConv, double Tconv, double cVconv):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mCVconv(cVconv) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    T = y/mTconv;
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
         const QuadraticInterpolator& epsMinInterp,
         const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mVelConv(velConv),
    mEpsMinInterp(epsMinInterp),
    mTinterp(Tinterp) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    const auto eps0 = mEpsMinInterp(x);
    T = mTinterp(x, y - eps0)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return cs * mVelConv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mVelConv;
  const QuadraticInterpolator& mEpsMinInterp;
  const BiQuadraticInterpolator& mTinterp;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// DPDrho(rho, eps)
//------------------------------------------------------------------------------
class Kfunc {
public:
  Kfunc(int matNum, double rhoConv, double Tconv, double Pconv,
        const QuadraticInterpolator& epsMinInterp,
        const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mPconv(Pconv),
    mEpsMinInterp(epsMinInterp),
    mTinterp(Tinterp) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    const auto eps0 = mEpsMinInterp(x);
    T = mTinterp(x, y - eps0)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return std::abs(rho * DPDR * mPconv);
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mPconv;
  const QuadraticInterpolator& mEpsMinInterp;
  const BiQuadraticInterpolator& mTinterp;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

//------------------------------------------------------------------------------
// s(rho, eps) (entropy)
//------------------------------------------------------------------------------
class sfunc {
public:
  sfunc(int matNum, double rhoConv, double Tconv, double Sconv,
        const QuadraticInterpolator& epsMinInterp,
        const BiQuadraticInterpolator& Tinterp):
    mMatNum(matNum),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mSconv(Sconv),
    mEpsMinInterp(epsMinInterp),
    mTinterp(Tinterp) {}
  double operator()(const double x, const double y) const {
    rho = x/mRhoConv;
    const auto eps0 = mEpsMinInterp(x);
    T = mTinterp(x, y - eps0)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return S * mSconv;
  }
private:
  mutable int mMatNum;
  double mRhoConv, mTconv, mSconv;
  const QuadraticInterpolator& mEpsMinInterp;
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
      const MaterialPressureMinType minPressureType,
      const bool useInterpolation):
  SolidEquationOfState<Dimension>(get_aneos_referencedensity_(const_cast<int*>(&materialNumber)),  // not in the right units yet!
                                  0.0,                                           // dummy etamin
                                  std::numeric_limits<double>::max(),            // dummy etamax
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minPressureType),
  mUseInterpolation(useInterpolation),
  mMaterialNumber(materialNumber),
  mNumRhoVals(numRhoVals),
  mNumTvals(numTvals),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax),
  mTmin(Tmin),
  mTmax(Tmax),
  mEpsMin(std::numeric_limits<double>::min()),
  mEpsMax(std::numeric_limits<double>::max()),
  mExternalPressure(externalPressure),
  mEpsMinInterp(),
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

  // Build the interpolator for the mininum specific energy along the min density line.
  // Also find the maximum specific energy.
  const auto Feps = epsFunc(mMaterialNumber, mRhoConv, mTconv, mEconv);
  const auto drho = (mRhoMax - mRhoMin)/(mNumRhoVals - 1u);
  CHECK(drho > 0.0);
  // [rhoMin:rhoMax], Tmin
  vector<double> epsMinVals(mNumRhoVals);
  for (auto i = 0u; i < mNumRhoVals; ++i) {
    epsMinVals[i] = Feps(mRhoMin + i*drho, mTmin);
    const auto epsi = Feps(mRhoMin + i*drho, mTmax);
    mEpsMax = min(mEpsMax, epsi);
  }
  mEpsMinInterp.initialize(mTmin, mTmax, epsMinVals);

  // Build the biquadratic interpolation function for eps(rho, T)
  auto t0 = clock();
  mEpsInterp.initialize(mRhoMin, mRhoMax,
                        mTmin, mTmax,
                        mNumRhoVals, mNumTvals, Feps);
  if (Process::getRank() == 0) cout << "ANEOS: Time to build epsInterp: " << double(clock() - t0)/CLOCKS_PER_SEC << endl;

  // Now the hard inversion method for looking up T(rho, eps)
  t0 = clock();
  const auto Ftemp = Tfunc(mTmin, mTmax, mEpsMinInterp, mEpsInterp, Feps);
  mTinterp.initialize(mRhoMin, mRhoMax,
                      0.0, mEpsMax,
                      mNumRhoVals, mNumTvals, Ftemp);
  if (Process::getRank() == 0) cout << "ANEOS: Time to build Tinterp: " << double(clock() - t0)/CLOCKS_PER_SEC << endl;

  // And finally the interpolators for most of our derived quantities
  t0 = clock();
  const auto Fpres = Pfunc(mMaterialNumber, mRhoConv, mTconv, mPconv, mEpsMinInterp, mTinterp);
  mPinterp.initialize(mRhoMin, mRhoMax,
                      0.0, mEpsMax,
                      mNumRhoVals, mNumTvals, Fpres);
  if (Process::getRank() == 0) cout << "ANEOS: Time to build Pinterp: " << double(clock() - t0)/CLOCKS_PER_SEC << endl;

  t0 = clock();
  const auto Fcv = cVfunc(mMaterialNumber, mRhoConv, mTconv, mCVconv);
  mCVinterp.initialize(mRhoMin, mRhoMax,
                       mTmin, mTmax,
                       mNumRhoVals, mNumTvals, Fcv);
  if (Process::getRank() == 0) cout << "ANEOS: Time to build CVinterp: " << double(clock() - t0)/CLOCKS_PER_SEC << endl;

  t0 = clock();
  const auto Fcs = csfunc(mMaterialNumber, mRhoConv, mTconv, mVelConv, mEpsMinInterp, mTinterp);
  mCSinterp.initialize(mRhoMin, mRhoMax,
                       0.0, mEpsMax,
                       mNumRhoVals, mNumTvals, Fcs);
  if (Process::getRank() == 0) cout << "ANEOS: Time to build CSinterp: " << double(clock() - t0)/CLOCKS_PER_SEC << endl;

  t0 = clock();
  const auto FK = Kfunc(mMaterialNumber, mRhoConv, mTconv, mPconv, mEpsMinInterp, mTinterp);
  mKinterp.initialize(mRhoMin, mRhoMax,
                      0.0, mEpsMax,
                      mNumRhoVals, mNumTvals, FK);
  if (Process::getRank() == 0) cout << "ANEOS: Time to build Kinterp: " << double(clock() - t0)/CLOCKS_PER_SEC << endl;

  t0 = clock();
  const auto Fs = sfunc(mMaterialNumber, mRhoConv, mTconv, mSconv, mEpsMinInterp, mTinterp);
  mSinterp.initialize(mRhoMin, mRhoMax,
                      0.0, mEpsMax,
                      mNumRhoVals, mNumTvals, Fs);
  if (Process::getRank() == 0) cout << "ANEOS: Time to build Sinterp: " << double(clock() - t0)/CLOCKS_PER_SEC << endl;
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
  double P;
  if (mUseInterpolation) {

    P = mPinterp(massDensity, specificThermalEnergy - mEpsMinInterp(massDensity));

  } else {

    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    P *= mPconv;

  }
  return this->applyPressureLimits(P - mExternalPressure);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  const auto eps0 = mEpsMinInterp(massDensity);
  const auto eps = specificThermalEnergy - eps0;
  if (not mUseInterpolation or (eps < 0.0 or eps > mEpsMax)) {
    const auto Feps = epsFunc(mMaterialNumber, mRhoConv, mTconv, mEconv);
    const auto Ftemp = Tfunc(mTmin, mTmax, mEpsMinInterp, mEpsInterp, Feps, true);
    return Ftemp(massDensity, specificThermalEnergy);
  } else {
    return mTinterp(massDensity, eps);
  }
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  if (mUseInterpolation) {
    return mEpsInterp(massDensity, temperature);
  } else {
    auto T = temperature/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return eps * mEconv;
  }    
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  if (mUseInterpolation) {
    return mCVinterp(massDensity, temperature);
  } else {
    auto T = temperature/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return cV * mCVconv;
  }    
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  // If we're out of bounds, extrapolate as though a gamma-law gas
  if (mUseInterpolation) {
    return mCSinterp(massDensity, specificThermalEnergy - mEpsMinInterp(massDensity));
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return cs * mVelConv;
  }    
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  const auto Ti = mTinterp(massDensity, specificThermalEnergy);     // this->temperature(massDensity, specificThermalEnergy);
  const auto cvi = mCVinterp(massDensity, Ti);                      // this->specificHeat(massDensity, Ti);
  const auto nDen = massDensity/mAtomicWeight;
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
  if (mUseInterpolation) {
    return mKinterp(massDensity, specificThermalEnergy - mEpsMinInterp(massDensity));
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return  std::abs(rho*DPDR) * mPconv;
  }
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ANEOS<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  if (mUseInterpolation) {
    return mSinterp(massDensity, specificThermalEnergy - mEpsMinInterp(massDensity));
  } else {
    auto T = this->temperature(massDensity, specificThermalEnergy)/mTconv;
    auto rho = massDensity/mRhoConv;
    double P, eps, S, cV, DPDT, DPDR, cs;
    call_aneos_(const_cast<int*>(&mMaterialNumber), &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    return S * mSconv;
  }    
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
double
ANEOS<Dimension>::
epsMin() const {
  return mEpsMin;
}

template<typename Dimension>
double
ANEOS<Dimension>::
epsMax() const {
  return mEpsMax;
}

template<typename Dimension>
bool
ANEOS<Dimension>::
useInterpolation() const {
  return mUseInterpolation;
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

