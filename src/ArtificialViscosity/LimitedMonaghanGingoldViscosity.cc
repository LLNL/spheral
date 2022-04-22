//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#include "LimitedMonaghanGingoldViscosity.hh"
#include "Boundary/Boundary.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// limiter for velocity projection.
//------------------------------------------------------------------------------
//double limiterBJ(const double x) {
//  if (x > 0.0) {
//    return min(1.0, 4.0/(x + 1.0)*min(1.0, x));  // Barth-Jesperson
//  } else {
//    return 0.0;
//  }
//}
//
//double limiterMC(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*min(2.0*x, min(0.5*(1.0 + x), 2.0));   // monotonized central
//  } else {
//    return 0.0;
//  }
//}

double limiterVL(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*2.0*x/(1.0 + x);                       // van Leer
  } else {
    return 0.0;
  }
}

//double limiterMM(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*min(1.0, x);                           // minmod
//  } else {
//    return 0.0;
//  }
//}
//
//double limiterSB(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*max(min(2.0*x, 1.0), min(x, 2.0));    // superbee
//  } else {
//    return 0.0;
//  }
//}

// template<typename Vector, typename Tensor>
// double limiterConservative(const Vector& vi, const Vector& vj,
//                            const Vector& xi, const Vector& xj,
//                            const Tensor& DvDxi, const Tensor& DvDxj) {
//   const auto xji = xj - xi;
//   const auto vji = vj - vi;
//   const auto di = DvDxi.dot(xji);
//   const auto dj = DvDxj.dot(xji);
//   if (di.dot(dj) <= 0.0 or
//       di.dot(vji) <= 0.0 or
//       dj.dot(vji) <= 0.0) return 0.0;
//   const auto vjimag = vji.magnitude();
//   const auto dimag = di.magnitude();
//   const auto djmag = dj.magnitude();
//   return min(1.0, min(abs(vjimag*safeInv(dimag)), abs(vjimag*safeInv(djmag))));
// }

//double limiterConservative(const double vji, const double deltavi, const double deltavj) {
//  if (deltavi*deltavj <= 0.0 or
//      deltavi*vji <= 0.0 or
//      deltavj*vji <= 0.0) return 0.0;
//  return min(1.0, min(abs(vji*safeInv(deltavi)), abs(vji*safeInv(deltavj))));
//}

}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
LimitedMonaghanGingoldViscosity<Dimension>::
LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                const Scalar Cquadratic,
                                const bool linearInExpansion,
                                const bool quadraticInExpansion,
                                const Scalar etaCritFrac,
                                const Scalar etaFoldFrac):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, 
                                      linearInExpansion, quadraticInExpansion),
  mEtaCritFrac(etaCritFrac),
  mEtaFoldFrac(etaFoldFrac),
  mEtaCrit(0.0),
  mEtaFold(1.0),
  mGradVel(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
LimitedMonaghanGingoldViscosity<Dimension>::
~LimitedMonaghanGingoldViscosity() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LimitedMonaghanGingoldViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           const State<Dimension>& state,
           const StateDerivatives<Dimension>& derivs,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
           const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  // Let the base class do it's thing.
  ArtificialViscosity<Dimension>::initialize(dataBase, state, derivs, boundaryBegin, boundaryEnd, time, dt, W);

  // Cache the last velocity gradient for use during the step.
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  mGradVel = DvDx;
  mGradVel.copyFields();

  // // If any points are flagged as surface, force zero velocity gradient.
  // if (state.fieldNameRegistered(HydroFieldNames::surfacePoint)) {
  //   const auto surface = state.fields(HydroFieldNames::surfacePoint, 0);
  //   // const auto m0 = state.fields(HydroFieldNames::m0_Limited, 0.0);
  //   const auto numNodeLists = mGradVel.size();
  //   CHECK(surfacePoint.size() == numNodeLists);
  //   for (auto k = 0; k < numNodeLists; ++k) {
  //     const auto nk = mGradVel[k]->numInternalElements();
  //     for (auto i = 0; i < nk; ++i) {
  //       if (surface(k,i) != 0) mGradVel(k,i).Zero();
  //       // const auto m0i = min(m0(k,i), 1.0/m0(k,i));
  //       // mGradVel(k,i) *= std::max(0.0, 2.0*m0i - 1.0);
  //     }
  //   }
  // }

  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr < boundaryEnd;
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mGradVel);
  }
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Store the eta_crit value based on teh nodes perh smoothing scale.
  const double nPerh = dynamic_cast<const FluidNodeList<Dimension>&>(mGradVel[0]->nodeList()).nodesPerSmoothingScale();
  mEtaCrit = mEtaCritFrac/nPerh;
  mEtaFold = mEtaFoldFrac/nPerh;
  CHECK(mEtaFold > 0.0);
}

//------------------------------------------------------------------------------
// The required method to compute the artificial viscous P/rho^2.
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
     typename Dimension::Tensor>
LimitedMonaghanGingoldViscosity<Dimension>::
Piij(const unsigned nodeListi, const unsigned i, 
     const unsigned nodeListj, const unsigned j,
     const Vector& xi,
     const Vector& etai,
     const Vector& vi,
     const Scalar rhoi,
     const Scalar csi,
     const SymTensor& Hi,
     const Vector& xj,
     const Vector& etaj,
     const Vector& vj,
     const Scalar rhoj,
     const Scalar csj,
     const SymTensor& Hj) const {

  double Cl = this->mClinear;
  double Cq = this->mCquadratic;
  const double eps2 = this->mEpsilon2;
  const bool linearInExp = this->linearInExpansion();
  const bool quadInExp = this->quadraticInExpansion();
  const bool balsaraShearCorrection = this->mBalsaraShearCorrection;
  const Tensor& DvDxi = mGradVel(nodeListi, i);
  const Tensor& DvDxj = mGradVel(nodeListj, j);

  // Grab the FieldLists scaling the coefficients.
  // These incorporate things like the Balsara shearing switch or Morris & Monaghan time evolved
  // coefficients.
  const Scalar fCli = this->mClMultiplier(nodeListi, i);
  const Scalar fCqi = this->mCqMultiplier(nodeListi, i);
  const Scalar fClj = this->mClMultiplier(nodeListj, j);
  const Scalar fCqj = this->mCqMultiplier(nodeListj, j);
  Cl *= 0.5*(fCli + fClj);
  Cq *= 0.5*(fCqi + fCqj);

  // Are we applying the shear corrections?
  Scalar fshear = 1.0;
  if (balsaraShearCorrection) {
    const Scalar csneg = this->negligibleSoundSpeed();
    const Scalar hiinv = Hi.Trace()/Dimension::nDim;
    const Scalar hjinv = Hj.Trace()/Dimension::nDim;
    const Scalar ci = max(csneg, csi);
    const Scalar cj = max(csneg, csj);
    const Scalar fi = this->curlVelocityMagnitude(DvDxi)/(this->curlVelocityMagnitude(DvDxi) + abs(DvDxi.Trace()) + eps2*ci*hiinv);
    const Scalar fj = this->curlVelocityMagnitude(DvDxj)/(this->curlVelocityMagnitude(DvDxj) + abs(DvDxj.Trace()) + eps2*cj*hjinv);
    fshear = min(fi, fj);
  }
/*
  else{
    const Tensor Shi = 0.5*(DvDxi+DvDxi.Transpose())-(1.0/Dimension::nDim)*Tensor::one*DvDxi.Trace();
    const Tensor Shj = 0.5*(DvDxj+DvDxj.Transpose())-(1.0/Dimension::nDim)*Tensor::one*DvDxj.Trace();
    const Scalar csneg = this->negligibleSoundSpeed();
    const Scalar hiinv = Hi.Trace()/Dimension::nDim;
    const Scalar hjinv = Hj.Trace()/Dimension::nDim;
    const Scalar ci = max(csneg, csi);
    const Scalar cj = max(csneg, csj);
    //const Tensor DivVi = (1.0/Dimension::nDim)*Tensor::one*DvDxi.Trace();
    //const Tensor DivVj = (1.0/Dimension::nDim)*Tensor::one*DvDxj.Trace();
    const Scalar DivVi = DvDxi.Trace();
    const Scalar DivVj = DvDxj.Trace();
    //const Scalar fi = (DivVi*(DivVi.Transpose())).Trace()/((DvDxi*(DvDxi.Transpose())).Trace() + eps2*ci*hiinv);
    //const Scalar fj = (DivVj*(DivVj.Transpose())).Trace()/((DvDxj*(DvDxj.Transpose())).Trace() + eps2*cj*hjinv);
    //const Scalar fi = (DivVi*(DivVi.Transpose())).Trace()/((DvDxi*(DvDxi.Transpose())).Trace() + 1e-30);
    //const Scalar fj = (DivVj*(DivVj.Transpose())).Trace()/((DvDxj*(DvDxj.Transpose())).Trace() + 1e-30);
    const Scalar fi = DivVi*DivVi*safeInv((Shi*Shi.Transpose()).Trace()+DivVi*DivVi);
    const Scalar fj = DivVj*DivVj*safeInv((Shj*Shj.Transpose()).Trace()+DivVj*DivVj);
    fshear = min(fi, fj);
    fshear = 1.0;
  }
*/

  // Compute the corrected velocity difference.
  // Vector vij = vi - vj;
  const auto xij = 0.5*(xi - xj);
  const auto gradi = (DvDxi.dot(xij)).dot(xij);
  const auto gradj = (DvDxj.dot(xij)).dot(xij);
  const auto ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  const auto rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  CHECK(min(ri, rj) <= 1.0);
  // const Scalar phi = limiterMM(min(ri, rj));
  auto phi = limiterVL(min(ri, rj));
  
  // const auto xjihat = -xij.unitVector();
  // auto phi = limiterConservative((vj - vi).dot(xjihat), (DvDxi*xjihat).dot(xjihat), (DvDxj*xjihat).dot(xjihat));

  // If the points are getting too close, we let the Q come back full force.
  const auto etaij = min(etai.magnitude(), etaj.magnitude());
  // phi *= (etaij2 < etaCrit2 ? 0.0 : 1.0);
  // phi *= min(1.0, etaij2*etaij2/(etaCrit2etaCrit2));
  if (etaij < mEtaCrit) {
    phi *= exp(-FastMath::square((etaij - mEtaCrit)/mEtaFold));
  }

  // "Mike" method.
  const Vector vi1 = vi - phi*DvDxi*xij;
  const Vector vj1 = vj + phi*DvDxj*xij;

  // const Vector vi1 = vi - phi*DvDxi*xij;
  // const Vector vj1 = vj + phi*DvDxj*xij;
  
  const Vector vij = vi1 - vj1;
  
  // Compute mu.
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

  // The artificial internal energy.
  const Scalar ei = fshear*(-Cl*csi*(linearInExp    ? mui                : min(0.0, mui)) +
                             Cq    *(quadInExp      ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui))));
  const Scalar ej = fshear*(-Cl*csj*(linearInExp    ? muj                : min(0.0, muj)) +
                             Cq    *(quadInExp      ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj))));
  CHECK2(ei >= 0.0 or (linearInExp or quadInExp), ei << " " << csi << " " << mui);
  CHECK2(ej >= 0.0 or (linearInExp or quadInExp), ej << " " << csj << " " << muj);

  // Now compute the symmetrized artificial viscous pressure.
  return make_pair(ei/rhoi*Tensor::one,
                   ej/rhoj*Tensor::one);
}

//------------------------------------------------------------------------------
// etaCritFrac
//------------------------------------------------------------------------------
template<typename Dimension>
double
LimitedMonaghanGingoldViscosity<Dimension>::
etaCritFrac() const {
  return mEtaCritFrac;
}

template<typename Dimension>
void
LimitedMonaghanGingoldViscosity<Dimension>::
etaCritFrac(double val) {
  VERIFY(val >= 0.0);
  mEtaCritFrac = val;
}

//------------------------------------------------------------------------------
// etaFoldFrac
//------------------------------------------------------------------------------
template<typename Dimension>
double
LimitedMonaghanGingoldViscosity<Dimension>::
etaFoldFrac() const {
  return mEtaFoldFrac;
}

template<typename Dimension>
void
LimitedMonaghanGingoldViscosity<Dimension>::
etaFoldFrac(double val) {
  VERIFY(val > 0.0);
  mEtaFoldFrac = val;
}

}
