//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "MonaghanGingoldViscosity.hh"
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

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldViscosity<Dimension>::
MonaghanGingoldViscosity(const Scalar Clinear,
                         const Scalar Cquadratic,
                         const bool linearInExpansion,
                         const bool quadraticInExpansion):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic),
  mLinearInExpansion(linearInExpansion),
  mQuadraticInExpansion(quadraticInExpansion) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldViscosity<Dimension>::
~MonaghanGingoldViscosity() {
}


//------------------------------------------------------------------------------
// The required method to compute the artificial viscous P/rho^2.
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
     typename Dimension::Tensor>
MonaghanGingoldViscosity<Dimension>::
Piij(const unsigned nodeListi, const unsigned i, 
     const unsigned nodeListj, const unsigned j,
     const Vector& /*xi*/,
     const Vector& etai,
     const Vector& vi,
     const Scalar rhoi,
     const Scalar csi,
     const SymTensor& /*Hi*/,
     const Vector& /*xj*/,
     const Vector& etaj,
     const Vector& vj,
     const Scalar rhoj,
     const Scalar csj,
     const SymTensor& /*Hj*/) const {

  double Cl = this->mClinear;
  double Cq = this->mCquadratic;
  const double eps2 = this->mEpsilon2;

  // Grab the FieldLists scaling the coefficients.
  // These incorporate things like the Balsara shearing switch or Morris & Monaghan time evolved
  // coefficients.
  const auto fCli = this->mClMultiplier(nodeListi, i);
  const auto fCqi = this->mCqMultiplier(nodeListi, i);
  const auto fClj = this->mClMultiplier(nodeListj, j);
  const auto fCqj = this->mCqMultiplier(nodeListj, j);
  const auto fshear = 0.5*(this->mShearCorrection(nodeListi, i) + this->mShearCorrection(nodeListj, j));
  Cl *= 0.5*(fCli + fClj)*fshear;
  Cq *= 0.5*(fCqi + fCqj)*fshear;

  // Scalar fshear = 1.0;
  // Scalar fsheari = fshear;
  // Scalar fshearj = fshear;
  // const Tensor& DvDxi = mGradVel(nodeListi, i);
  // const Tensor& DvDxj = mGradVel(nodeListj, j);
  // if (balsaraShearCorrection) {
  //   const Scalar csneg = this->negligibleSoundSpeed();
  //   const Scalar hiinv = Hi.Trace()/Dimension::nDim;
  //   const Scalar hjinv = Hj.Trace()/Dimension::nDim;
  //   const Scalar ci = max(csneg, csi);
  //   const Scalar cj = max(csneg, csj);
  //   const Scalar fi = abs(DvDxi.Trace())/(this->curlVelocityMagnitude(DvDxi) + abs(DvDxi.Trace()) + eps2*ci*hiinv);
  //   const Scalar fj = abs(DvDxj.Trace())/(this->curlVelocityMagnitude(DvDxj) + abs(DvDxj.Trace()) + eps2*cj*hjinv);
  //   fshear = min(fi, fj);
  //   //fsheari = fi;
  //   //fshearj = fj;
  //   fsheari = fshear;
  //   fshearj = fshear;
  // }

  // Compute mu.
  const auto vij = vi - vj;
  const auto mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const auto muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

  // The artificial internal energy.
  const auto ei = -Cl*csi*(mLinearInExpansion    ? mui                : min(0.0, mui)) +
                   Cq    *(mQuadraticInExpansion ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui)));
  const auto ej = -Cl*csj*(mLinearInExpansion    ? muj                : min(0.0, muj)) +
                   Cq    *(mQuadraticInExpansion ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj)));
  CHECK2(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ei << " " << csi << " " << mui);
  CHECK2(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ej << " " << csj << " " << muj);

  // Now compute the symmetrized artificial viscous pressure.
  return make_pair(ei/rhoi*Tensor::one,
                   ej/rhoj*Tensor::one);
}

//------------------------------------------------------------------------------
// linearInExpansion
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MonaghanGingoldViscosity<Dimension>::
linearInExpansion() const {
  return mLinearInExpansion;
}

template<typename Dimension>
void
MonaghanGingoldViscosity<Dimension>::
linearInExpansion(const bool x) {
  mLinearInExpansion = x;
}

//------------------------------------------------------------------------------
// quadraticInExpansion
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MonaghanGingoldViscosity<Dimension>::
quadraticInExpansion() const {
  return mQuadraticInExpansion;
}

template<typename Dimension>
void
MonaghanGingoldViscosity<Dimension>::
quadraticInExpansion(const bool x) {
  mQuadraticInExpansion = x;
}

}
