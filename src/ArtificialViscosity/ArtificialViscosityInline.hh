#include "DataBase/DataBase.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the linear coefficient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ArtificialViscosity<Dimension>::
Cl() const {
  return mClinear;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
Cl(typename Dimension::Scalar Cl) {
  mClinear = Cl;
}

//------------------------------------------------------------------------------
// Access the quadratic coefficient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ArtificialViscosity<Dimension>::
Cq() const {
  return mCquadratic;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
Cq(typename Dimension::Scalar Cq) {
  mCquadratic = Cq;
}


//------------------------------------------------------------------------------
// Access QcorrectionOrder
//------------------------------------------------------------------------------
template<typename Dimension>
inline
RKOrder
ArtificialViscosity<Dimension>::QcorrectionOrder() const {
  return mQcorrectionOrder;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
QcorrectionOrder(const RKOrder order) {
  mQcorrectionOrder = order;
}


//------------------------------------------------------------------------------
// Access sigma.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
ArtificialViscosity<Dimension>::
sigma() const {
  return mSigma;
}

//------------------------------------------------------------------------------
// Toggle for the Balsara shear correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ArtificialViscosity<Dimension>::
balsaraShearCorrection() const {
  return mBalsaraShearCorrection;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
balsaraShearCorrection(bool value) {
  mBalsaraShearCorrection = value;
}

//------------------------------------------------------------------------------
// Access the FieldLists of multipliers.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::Scalar>&
ArtificialViscosity<Dimension>::
ClMultiplier() {
  return mClMultiplier;
}
    
template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::Scalar>&
ArtificialViscosity<Dimension>::
CqMultiplier() {
  return mCqMultiplier;
}
    
template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::Scalar>&
ArtificialViscosity<Dimension>::
shearCorrection() {
  return mShearCorrection;
}
    
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
ArtificialViscosity<Dimension>::
ClMultiplier() const {
  return mClMultiplier;
}
    
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
ArtificialViscosity<Dimension>::
CqMultiplier() const {
  return mCqMultiplier;
}
    
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
ArtificialViscosity<Dimension>::
shearCorrection() const {
  return mShearCorrection;
}
    
//------------------------------------------------------------------------------
// Calculate the curl of the velocity given the stress tensor.
//------------------------------------------------------------------------------
template<>
inline
Dim<1>::Scalar
ArtificialViscosity< Dim<1> >::
curlVelocityMagnitude(const Dim<1>::Tensor& /*DvDx*/) const {
  return 0.0;
}

template<>
inline
Dim<2>::Scalar
ArtificialViscosity< Dim<2> >::
curlVelocityMagnitude(const Dim<2>::Tensor& DvDx) const {
  return std::abs(DvDx.yx() - DvDx.xy());
}

template<>
inline
Dim<3>::Scalar
ArtificialViscosity< Dim<3> >::
curlVelocityMagnitude(const Dim<3>::Tensor& DvDx) const {
  return sqrt(FastMath::square(DvDx.zy() - DvDx.yz()) +
	      FastMath::square(DvDx.xz() - DvDx.zx()) +
	      FastMath::square(DvDx.yx() - DvDx.xy()));
}

//------------------------------------------------------------------------------
// Access the gradient of the velocity divergence.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
ArtificialViscosity<Dimension>::
gradDivVelocity() const {
  return mGradDivVelocity;
}

//------------------------------------------------------------------------------
// Access the switch to activate the limiter.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ArtificialViscosity<Dimension>::
limiter() const {
  return mLimiterSwitch;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
limiter(bool value) {
  mLimiterSwitch = value;
  if (mLimiterSwitch == true) calculateSigma(true);
}

//------------------------------------------------------------------------------
// Access the neglibile sound speed parameter for the limiter.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ArtificialViscosity<Dimension>::
negligibleSoundSpeed() const {
  return mNegligibleSoundSpeed;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
negligibleSoundSpeed(typename Dimension::Scalar val) {
  REQUIRE(val > 0.0);
  mNegligibleSoundSpeed = val;
}

//------------------------------------------------------------------------------
// Access the sound speed multiplier for the limiter.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ArtificialViscosity<Dimension>::
csMultiplier() const {
  return mCsMultiplier;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
csMultiplier(typename Dimension::Scalar val) {
  REQUIRE(val > 0.0);
  mCsMultiplier = val;
}

//------------------------------------------------------------------------------
// Access the energy multiplier for the limiter.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ArtificialViscosity<Dimension>::
energyMultiplier() const {
  return mEnergyMultiplier;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
energyMultiplier(typename Dimension::Scalar val) {
  REQUIRE(val > 0.0);
  mEnergyMultiplier = val;
}

//------------------------------------------------------------------------------
// Access the epsilon2 term.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ArtificialViscosity<Dimension>::
epsilon2() const {
  return mEpsilon2;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
epsilon2(typename Dimension::Scalar epsilon2) {
  mEpsilon2 = epsilon2;
}

//------------------------------------------------------------------------------
// Compute the limiter function for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Tensor
ArtificialViscosity<Dimension>::
calculateLimiter(const Vector& /*vi*/,
                 const Vector& /*vj*/,
                 const Scalar  ci,
                 const Scalar  /*cj*/,
                 const Scalar  hi,
                 const Scalar  /*hj*/,
                 const int nodeListID,
                 const int nodeID) const {

  using std::max;
  using std::min;
  using std::abs;

  REQUIRE(negligibleSoundSpeed() > 0.0);
  REQUIRE(epsilon2() > 0.0);

//------------------------------------------------------------------------------
  const Vector shati = shockDirection(ci, hi, nodeListID, nodeID);
  return shati.dyad(shati);

// //------------------------------------------------------------------------------
//   const Vector shati = shockDirection(ci, hi, nodeI);
//   const Vector shatj = shockDirection(cj, hj, nodeJ);
//   const Scalar fi = max(0.0, -(shati.dot(shatj)));
//   return max(mEnergyMultiplier, fi*fi);

// //------------------------------------------------------------------------------
//   const Vector shati = shockDirection(ci, hi, nodeI);
//   const Vector shatj = shockDirection(cj, hj, nodeJ);
//   const Scalar fi = shati.dot(shatj);
//   return fi*fi;

// //------------------------------------------------------------------------------
//   const Vector shati = shockDirection(ci, hi, nodeI);
//   const Vector shatj = shockDirection(cj, hj, nodeJ);
//   const Scalar mag = max(0.0, -(shati.dot(shatj)));
//   return mag*mIdentity;

// //------------------------------------------------------------------------------
//   const Vector shati = shockDirection(ci, hi, nodeI);
//   const Vector shatj = shockDirection(cj, hj, nodeJ);
// //   const Scalar mag = max(0.0, -(shati.dot(shatj)));
//   const Scalar mag = shati.dot(shatj);
//   const Vector shocki = abs(mag)*shati;
//   const Vector shockj = abs(mag)*shatj;
//   if (mag < 0.0) {
//     return 0.5*(shocki.dyad(shocki) - shockj.dyad(shockj));
//   } else {
//     return 0.5*(shocki.dyad(shocki) + shockj.dyad(shockj));
//   }

// //------------------------------------------------------------------------------
//   const Scalar mi = mMass(nodeI);
//   const Scalar mj = mMass(nodeJ);
//   const Scalar csi = max(negligibleSoundSpeed(), epsilon2()*ci);
//   CHECK(csi > 0.0);

//   Vector gDivVi = mGradDivVelocity(nodeI);
//   const Vector gDivVhati = gDivVi.unitVector();
//   const Scalar divVi = abs(mSigma(nodeI).Trace());
//   const Scalar vis = vi.dot(gDivVhati);
//   const Scalar vjs = vj.dot(gDivVhati);
//   const Scalar dKEmax = 0.5*abs(mi*vis*vis - mj*vjs*vjs);
//   const Scalar dKEq = mi*hi*hi*(2.0*Cq()*hi*divVi + Cl()*ci + csi);
//   CHECK(dKEq > 0.0);
//   Vector sHat = gDivVhati*min(1.0, mEnergyMultiplier*dKEmax/dKEq);

//   const Scalar weight = gDivVhati.magnitude2();
//   CHECK(fuzzyGreaterThanOrEqual(weight, 0.0) &&
//         fuzzyLessThanOrEqual(weight, 1.0));

//   return sHat.dyad(sHat) + (1.0 - weight)*mShearMultiplier(nodeI)*mIdentity;

// //------------------------------------------------------------------------------
//   const Scalar csi = max(negligibleSoundSpeed(), epsilon2()*ci);
//   const Scalar csOverhi2 = csi/(hi*hi);
//   CHECK(csOverhi2 > 0.0);

//   Vector gDivVi = mGradDivVelocity(nodeI);
//   const Vector gDivVhati = gDivVi.unitVector();
//   const Scalar divVi = abs(mSigma(nodeI).Trace());
//   const Scalar thpt = vij.dot(gDivVhati) + csi;
//   const Scalar KE = 0.5*thpt*thpt;
//   const Scalar ack = hi*hi*(2.0*divVi*hi + csi);
//   CHECK(ack > 0.0);
//   Vector sHat = gDivVi/max(gDivVi.magnitude() + 1.0e-10,
//                            mEnergyMultiplier*KE/ack);

// //   const Scalar weight = gDivVhati.magnitude2();
// //   CHECK(fuzzyGreaterThanOrEqual(weight, 0.0) &&
// //         fuzzyLessThanOrEqual(weight, 1.0));

//   return sHat.dyad(sHat); // + (1.0 - weight)*mShearMultiplier(nodeI)*mIdentity;

// //------------------------------------------------------------------------------
// //------------------------------------------------------------------------------
// //------------------------------------------------------------------------------
//   const Vector gDivVhati = shockDirection(rijUnit, ci, hi, nodeI);
//   const Scalar weight = gDivVhati.magnitude2();
//   CHECK(weight <= 1.0);
//   return gDivVhati.dyad(gDivVhati) + 
//     (1.0 - weight)*mShearMultiplier(nodeI)*mIdentity;

// //------------------------------------------------------------------------------
//   const Vector gDivVhati = shockDirection(rijUnit, ci, hi, nodeI);
//   const Scalar weight = gDivVhati.magnitude2();
//   CHECK(weight <= 1.0);
//   return mShearMultiplier(nodeI)*(gDivVhati.dyad(gDivVhati) + 
//                                   (1.0 - weight)*mIdentity);

}

//------------------------------------------------------------------------------
// Compute the shock direction indicator based on grad div velocity for the 
// given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
ArtificialViscosity<Dimension>::
shockDirection(const Scalar ci,
               const Scalar hi,
               const int nodeListID,
               const int nodeID) const {

  using std::max;
  using std::min;
  using std::abs;

  REQUIRE(hi > 0.0);
  REQUIRE(csMultiplier() > 0.0);
  REQUIRE(negligibleSoundSpeed() > 0.0);
  REQUIRE(epsilon2() > 0.0);
  REQUIRE(nodeListID >= 0 && nodeListID < mGradDivVelocity.size());
  REQUIRE(nodeID >= 0 && nodeID < mGradDivVelocity[nodeListID]->nodeListPtr()->numNodes());

  const Scalar csi = max(mNegligibleSoundSpeed, mCsMultiplier*ci);
  const Scalar csOverhi2 = csi/(hi*hi);
  CHECK(csOverhi2 > 0.0);

  Vector gDivVhati = (*mGradDivVelocity[nodeListID])(nodeID);
  const Scalar thpt = gDivVhati.magnitude() + csOverhi2;
  gDivVhati /= thpt;

//   const Scalar divVi = std::max(0.0, -(mSigma(nodeI).Trace()));
//   gDivVhati *= divVi/(divVi + csi/hi);

  ENSURE(fuzzyLessThanOrEqual(gDivVhati.magnitude2(), 1.0));
  return gDivVhati;
}

//------------------------------------------------------------------------------
// Access the switch to calculate the stress-strain tensor (sigma).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ArtificialViscosity<Dimension>::
calculateSigma() const {
  return mCalculateSigma;
}

template<typename Dimension>
inline
void
ArtificialViscosity<Dimension>::
calculateSigma(bool value) {
  mCalculateSigma = value;
}

//------------------------------------------------------------------------------
// Compute Del cross V based on the given stress-strain tensor.
//------------------------------------------------------------------------------
template<>
inline
Dim<1>::Scalar
ArtificialViscosity< Dim<1> >::
computeDelCrossVMagnitude(const Dim<1>::Tensor& /*sigma*/) const {
  return 0.0;
}

template<>
inline
Dim<2>::Scalar
ArtificialViscosity< Dim<2> >::
computeDelCrossVMagnitude(const Dim<2>::Tensor& sigma) const {
  return std::abs(sigma(1,0) - sigma(0,1));
}

template<>
inline
Dim<3>::Scalar
ArtificialViscosity< Dim<3> >::
computeDelCrossVMagnitude(const Dim<3>::Tensor& sigma) const {
  return sqrt(FastMath::square(sigma(2,1) - sigma(1,2)) +
              FastMath::square(sigma(2,0) - sigma(0,2)) +
              FastMath::square(sigma(1,0) - sigma(0,1)));
}

//------------------------------------------------------------------------------
// Helper to calculate the weighting for use in the sigma calculation.
//------------------------------------------------------------------------------
template<>
inline
Dim<1>::Vector
ArtificialViscosity<Dim<1> >::
sigmaWeighting(const Dim<1>::Vector&) const {
  return Dim<1>::Vector(1.0);
}

template<>
inline
Dim<2>::Vector
ArtificialViscosity<Dim<2> >::
sigmaWeighting(const Dim<2>::Vector& r) const {
  return Dim<2>::Vector(FastMath::square(r.x()),
                        FastMath::square(r.y()))/(r.magnitude2() + 1.0e-10);
}

template<>
inline
Dim<3>::Vector
ArtificialViscosity<Dim<3> >::
sigmaWeighting(const Dim<3>::Vector& r) const {
  return Dim<3>::Vector(FastMath::square(r.x()),
                        FastMath::square(r.y()),
                        FastMath::square(r.z()))/(r.magnitude2() + 1.0e-10);
}

}

