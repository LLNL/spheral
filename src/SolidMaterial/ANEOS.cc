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
  typedef typename boost::multi_array<double, 2> array_type;
  typedef typename array_type::array_view<1>::type slice_type;
  typedef typename array_type::const_array_view<1>::type const_slice_type;
  typedef boost::multi_array_types::index_range range;

public:
  Tfunc(const int matNum,
        double rhoMin, 
        double rhoMax, 
        double Tmin, 
        double Tmax,
        unsigned numRhoVals, 
        unsigned numTvals,
        double rhoConv,
        double Tconv, 
        double Econv,
        const BiQuadraticInterpolator& epsInterp):
    mNumRhoVals(numRhoVals),
    mNumTvals(numTvals),
    mRhoMin(rhoMin),
    mRhoMax(rhoMax),
    mTmin(Tmin),
    mTmax(Tmax),
    mRhoConv(rhoConv),
    mTconv(Tconv),
    mEconv(Econv),
    mdrho((mRhoMax - mRhoMin)/(mNumRhoVals - 1u)),
    mdT((mTmax - mTmin)/(mNumTvals - 1u)),
    mSTEvals(boost::extents[numRhoVals][numTvals]),
    mMatNum(matNum) {

    // Compute the tabulated eps(rho, T) values.
    double rhoi, Ti;
    for (auto i = 0; i < mNumRhoVals; ++i) {
      rhoi = min(mRhoMin + i*mdrho, mRhoMax);
      for (auto j = 0; j < mNumTvals; ++j) {
        Ti = min(mTmin + j*mdT, mTmax);
        mSTEvals[i][j] = epsInterp(rhoi, Ti);
      }
    }
  }

  double operator()(const Dim<2>::Vector& x) const {

    // Find the closest starting point on the tabulated surface.
    const auto rhoi = x[0], epsi = x[1];
    const auto irho0 = min(mNumRhoVals - 1, int(max(0.0, (rhoi - mRhoMin)/mdrho)));
    const_slice_type rho0_slice = mSTEvals[boost::indices[irho0][range(0, mNumTvals)]];
    const auto iT0 = max(0, min(int(mNumTvals - 1), bisectSearch(rho0_slice.begin(), rho0_slice.end(), epsi)));

    // Try to iterate the EOS until we're close to the desired specific thermal energy
    auto iter = 0u;
    auto epsgoal = epsi/mEconv;
    rho = rhoi/mRhoConv;
    T = (mTmin + iT0*mdT)/mTconv;
    call_aneos_(&mMatNum, &T, &rho,
                &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    while (abs(eps - epsgoal)/std::max(1e-15, abs(epsgoal)) > 1.0e-8 and ++iter < 20u) {
      T += 0.9*(epsgoal - eps)/cV;
      call_aneos_(&mMatNum, &T, &rho,
                  &P, &eps, &S, &cV, &DPDT, &DPDR, &cs);
    }

    return T * mTconv;
  }

private:
  int mNumRhoVals, mNumTvals;
  double mRhoMin, mRhoMax, mTmin, mTmax, mRhoConv, mTconv, mEconv, mdrho, mdT;
  array_type mSTEvals;
  mutable int mMatNum;
  mutable double rho, T, P, eps, S, cV, DPDT, DPDR, cs;
};

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
  const auto Ftemp = Tfunc(mMaterialNumber, mRhoMin, mRhoMax, mTmin, mTmax,
                           mNumRhoVals, mNumTvals, mRhoConv, mTconv, mEconv,
                           mEpsInterp);
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

