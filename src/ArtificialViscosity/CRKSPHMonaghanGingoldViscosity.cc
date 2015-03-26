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
// Sweby limiter
// \beta \in [1, 2], where \beta=2 corresponds to the superbee limit.
//------------------------------------------------------------------------------
double swebyLimiter(const double x, const double beta) {
  return max(0.0, max(min(beta*x, 1.0), min(x, beta)));
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
                               const bool quadraticInExpansion,
                               const Scalar beta):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, 
                                      linearInExpansion, quadraticInExpansion),
  mBeta(beta),
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
  const Scalar rj = gradj*safeInv(gradi);
  const Scalar ri = safeInv(rj);
  // const Scalar curli = this->curlVelocityMagnitude(DvDxi);
  // const Scalar curlj = this->curlVelocityMagnitude(DvDxj);
  // const Scalar divi = abs(DvDxi.Trace());
  // const Scalar divj = abs(DvDxj.Trace());
  // const Scalar betaij = min(2.0, 1.0 + 0.5*min(curli*safeInv(curli + divi), curlj*safeInv(curlj + divj)));
  const Scalar phii = swebyLimiter(ri, mBeta);
  const Scalar phij = swebyLimiter(rj, mBeta);

  // // const Scalar phii = max(0.0, min(2.0*ri, min(0.5*(1.0 + ri), 2.0))); // Van Leer (1)
  // // const Scalar phij = max(0.0, min(2.0*rj, min(0.5*(1.0 + rj), 2.0))); // Van Leer (1)
  // // const Scalar phii = max(0.0, min(1.0, (ri + abs(ri))/(1.0 + abs(ri)))); // Van Leer (2)
  // // const Scalar phij = max(0.0, min(1.0, (rj + abs(rj))/(1.0 + abs(rj)))); // Van Leer (2)
  // const Scalar phii = max(0.0, max(min(2.0*ri, 1.0), min(ri, 2.0))); // superbee
  // const Scalar phij = max(0.0, max(min(2.0*rj, 1.0), min(rj, 2.0))); // superbee
    // const Scalar phii = max(0.0, max(min(2.5*ri, 1.0), min(ri, 2.5)));  // Sweby
    // const Scalar phij = max(0.0, max(min(2.5*rj, 1.0), min(rj, 2.5)));  // Sweby
  // // const Scalar phii = max(0.0, min(1.0, 2.0*(ri + abs(ri))*safeInv(ri + 3.0))); // HQUICK
  // // const Scalar phij = max(0.0, min(1.0, 2.0*(rj + abs(rj))*safeInv(rj + 3.0))); // HQUICK
  // // const Scalar phii = max(0.0, min(1.0, ri)); // minmod
  // // const Scalar phij = max(0.0, min(1.0, rj)); // minmod
  // // const Scalar phii = (ri <= 0.0 ? 0.0 : ri*(3.0*ri + 1.0)*safeInv(FastMath::square(ri + 1.0))); // CHARM
  // // const Scalar phij = (rj <= 0.0 ? 0.0 : rj*(3.0*rj + 1.0)*safeInv(FastMath::square(rj + 1.0))); // CHARM
  // // const Scalar phii = max(0.0, min(ri, 2.0)); // Osher
  // // const Scalar phij = max(0.0, min(rj, 2.0)); // Osher
  
  // // "Mike" method.
  // const Vector vi1 = vi - phii*DvDxi*xij;
  // const Vector vj1 = vj + phij*DvDxj*xij;
  // vij = vi1 - vj1;
  
  // // "Nick" method.
  // const Vector vi1 = vi - DvDxi*xij;
  // const Vector vj1 = vj + DvDxj*xij;
  // vij = vij + min(phii,phij)*((vi1-vj1)-vij);
  
  // // "Nick" method Mark II. 
  // const Vector vi1 = vi - (DvDxi-(1.0/Dimension::nDim)*Tensor::one*DvDxi.Trace())*xij;
  // const Vector vj1 = vj + (DvDxj-(1.0/Dimension::nDim)*Tensor::one*DvDxj.Trace())*xij;

  // // "Nick" method Mark III. 
  // const Tensor DvDxii = DvDxi - (1.0/Dimension::nDim)*Tensor::one*DvDxi.Trace();
  // const Tensor DvDxjj = DvDxj - (1.0/Dimension::nDim)*Tensor::one*DvDxj.Trace();
  // const Scalar phi = min(phii, phij);
  // const Vector vi1 = vi - (phi*DvDxi + (1.0 - phi)*DvDxii)*xij;
  // const Vector vj1 = vj + (phi*DvDxj + (1.0 - phi)*DvDxjj)*xij;

  // "Mike" method Mark III. 
  const Scalar phi = min(phii, phij);
  const Vector vi1 = vi - phi*DvDxi*xij;
  const Vector vj1 = vj + phi*DvDxj*xij;

  // // "Mike" method Mark IV. 
  // const Scalar phiS = min(phiSi, phiSj);
  // const Scalar phiA = min(phiAi, phiAj);
  // const Vector vi1 = vi - (phiS*Si + phiA*Ai)*xij;
  // const Vector vj1 = vj + (phiS*Sj + phiA*Aj)*xij;

  //const Vector vi1 = vi - (DvDxi-(1.0/3.0)*Tensor::one*DvDxi.Trace())*xij*min(phii,phij);
  //const Vector vj1 = vj + (DvDxj-(1.0/3.0)*Tensor::one*DvDxj.Trace())*xij*min(phii,phij);

  //const Vector vi1 = vi - (0.5*(DvDxi+DvDxi.Transpose())-(1.0/3.0)*Tensor::one*DvDxi.Trace())*xij;
  //const Vector vj1 = vj + (0.5*(DvDxj+DvDxj.Transpose())-(1.0/3.0)*Tensor::one*DvDxj.Trace())*xij;

  const Vector vij1 = vi1 - vj1;

  //const Scalar phi = max(0.0, min(0.95, min(phii, phij)));
  //const Scalar phi = min(phii, phij);
  //vij = (1.0 - phi)*vij + phi*vij1;
  vij=vij1;
  
  // Compute mu.
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

  // The artificial internal energy.
  const Scalar ei = fshear*(-Cl*rvAlphaL(nodeListi,i)*csi*(linearInExp    ? mui                : min(0.0, mui)) +
                             Cq *rvAlphaQ(nodeListi,i)   *(quadInExp      ? -sgn(mui)*mui*mui : FastMath::square(min(0.0, mui)))) ;
  const Scalar ej = fshear*(-Cl*rvAlphaL(nodeListj,j)*csj*(linearInExp    ? muj                : min(0.0, muj)) +
                             Cq *rvAlphaQ(nodeListj,j)    *(quadInExp     ? -sgn(muj)*muj*muj : FastMath::square(min(0.0, muj))));
  CHECK2(ei >= 0.0 or (linearInExp or quadInExp), ei << " " << csi << " " << mui);
  CHECK2(ej >= 0.0 or (linearInExp or quadInExp), ej << " " << csj << " " << muj);

  // Now compute the symmetrized artificial viscous pressure.
  return make_pair(ei/rhoi*Tensor::one,
                   ej/rhoj*Tensor::one);
}

}
}
