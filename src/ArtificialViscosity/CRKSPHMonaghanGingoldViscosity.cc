//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "CRKSPHMonaghanGingoldViscosity.hh"
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
#include "CRKSPH/gradientCRKSPH.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using DataOutput::Restart;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::Neighbor;
using Material::EquationOfState;
using BoundarySpace::Boundary;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;

namespace {

//------------------------------------------------------------------------------
// limiter for velocity projection.
//------------------------------------------------------------------------------
double limiter(const double x) {
  if (x > 0.0) {
    // return min(1.0, 4.0/(x + 1.0)*min(1.0, x));  // Barth-Jesperson
    return 2.0/(1.0 + x)*
      // min(2.0*x, min(0.5*(1.0 + x), 2.0));   // monotonized central
      2.0*x/(1.0 + x);        // van Leer
      // min(1.0, x);            // minmod
  } else {
    return 0.0;
  }
}

}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHMonaghanGingoldViscosity<Dimension>::
CRKSPHMonaghanGingoldViscosity(const Scalar Clinear,
                               const Scalar Cquadratic,
                               const bool linearInExpansion,
                               const bool quadraticInExpansion):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, 
                                      linearInExpansion, quadraticInExpansion),
  mGradVel(FieldSpace::Reference) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHMonaghanGingoldViscosity<Dimension>::
~CRKSPHMonaghanGingoldViscosity() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHMonaghanGingoldViscosity<Dimension>::
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

  // Cache pointers to the velocity gradient.
  FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  mGradVel = DvDx;
}

//------------------------------------------------------------------------------
// The required method to compute the artificial viscous P/rho^2.
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
     typename Dimension::Tensor>
CRKSPHMonaghanGingoldViscosity<Dimension>::
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

  const double Cl = this->mClinear;
  const double Cq = this->mCquadratic;
  const double eps2 = this->mEpsilon2;
  const bool linearInExp = this->linearInExpansion();
  const bool quadInExp = this->quadraticInExpansion();
  const bool balsaraShearCorrection = this->mBalsaraShearCorrection;
  const FieldSpace::FieldList<Dimension, Scalar>& rvAlphaQ = this->reducingViscosityMultiplierQ();
  const FieldSpace::FieldList<Dimension, Scalar>& rvAlphaL = this->reducingViscosityMultiplierL();
  const Tensor& DvDxi = mGradVel(nodeListi, i);
  const Tensor& DvDxj = mGradVel(nodeListj, j);

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

  // Compute the corrected velocity difference.
  Vector vij = vi - vj;
  const Vector xij = 0.5*(xi - xj);
  // const SymTensor Si = DvDxi.Symmetric();
  // const SymTensor Sj = DvDxj.Symmetric();
  // const Tensor Ai = DvDxi.SkewSymmetric();
  // const Tensor Aj = DvDxj.SkewSymmetric();
  // const Scalar gradSi = (Si.dot(xij)).dot(xij);
  // const Scalar gradSj = (Sj.dot(xij)).dot(xij);
  // const Scalar gradAi = (Ai.dot(xij)).dot(xij);
  // const Scalar gradAj = (Aj.dot(xij)).dot(xij);
  // const Scalar rSi = gradSi*safeInv(gradSj);
  // const Scalar rSj = safeInv(rSi);
  // const Scalar rAi = gradAi*safeInv(gradAj);
  // const Scalar rAj = safeInv(rAi);
  // const Scalar phiSi = swebyLimiter(rSi, 1.5);
  // const Scalar phiSj = swebyLimiter(rSj, 1.5);
  // const Scalar phiAi = swebyLimiter(rAi, 2.0);
  // const Scalar phiAj = swebyLimiter(rAj, 2.0);

  const Scalar gradi = (DvDxi.dot(xij)).dot(xij);
  const Scalar gradj = (DvDxj.dot(xij)).dot(xij);
  const Scalar ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  const Scalar rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  // const Scalar curli = this->curlVelocityMagnitude(DvDxi);
  // const Scalar curlj = this->curlVelocityMagnitude(DvDxj);
  // const Scalar divi = abs(DvDxi.Trace());
  // const Scalar divj = abs(DvDxj.Trace());
  // const Scalar betaij = min(2.0, 1.0 + min(curli/max(1.0e-30, curli + divi), curlj/max(1.0e-30, curlj + divj)));

  // const Vector vij12 = 0.5*(vi + vj);
  // const Scalar phimax = min(1.0, abs(vij.dot(xij)*safeInv(vij12.dot(xij))));

  const Scalar phii = limiter(ri);
  const Scalar phij = limiter(rj);

  // "Mike" method.
  const Vector vi1 = vi - phii*DvDxi*xij;
  const Vector vj1 = vj + phij*DvDxj*xij;
  vij = vi1 - vj1;
  
  // Compute mu.
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

  // The artificial internal energy.
  const Scalar ei = fshear*(-Cl*rvAlphaL(nodeListi,i)*csi*(linearInExp    ? mui                : min(0.0, mui)) +
                             Cq *rvAlphaQ(nodeListi,i)   *(quadInExp      ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui)))) ;
  const Scalar ej = fshear*(-Cl*rvAlphaL(nodeListj,j)*csj*(linearInExp    ? muj                : min(0.0, muj)) +
                             Cq *rvAlphaQ(nodeListj,j)   *(quadInExp      ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj))));
  CHECK2(ei >= 0.0 or (linearInExp or quadInExp), ei << " " << csi << " " << mui);
  CHECK2(ej >= 0.0 or (linearInExp or quadInExp), ej << " " << csj << " " << muj);

  // Now compute the symmetrized artificial viscous pressure.
  return make_pair(ei/rhoi*Tensor::one,
                   ej/rhoj*Tensor::one);
}

}
}
