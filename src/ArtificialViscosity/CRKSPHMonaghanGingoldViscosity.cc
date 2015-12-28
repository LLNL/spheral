//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "CRKSPHMonaghanGingoldViscosity.hh"
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
#include "CRKSPH/computeCRKSPHMoments.hh"
#include "CRKSPH/computeCRKSPHCorrections.hh"
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
double limiterBJ(const double x) {
  if (x > 0.0) {
    return min(1.0, 4.0/(x + 1.0)*min(1.0, x));  // Barth-Jesperson
  } else {
    return 0.0;
  }
}

double limiterMC(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*min(2.0*x, min(0.5*(1.0 + x), 2.0));   // monotonized central
  } else {
    return 0.0;
  }
}

double limiterVL(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*2.0*x/(1.0 + x);                       // van Leer
  } else {
    return 0.0;
  }
}

double limiterMM(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*min(1.0, x);                           // minmod
  } else {
    return 0.0;
  }
}

double limiterSB(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*max(min(2.0*x, 1.0), min(x, 2.0));    // superbee
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
                               const bool quadraticInExpansion,
                               const Scalar etaCritFrac,
                               const Scalar etaFoldFrac):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, 
                                      linearInExpansion, quadraticInExpansion),
  mEtaCritFrac(etaCritFrac),
  mEtaFoldFrac(etaFoldFrac),
  mEtaCrit(0.0),
  mEtaFold(1.0),
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

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  // const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);

  const CRKSPHSpace::CRKOrder correctionOrder = this->QcorrectionOrder();

  // Make new correction fields cause do not want to mess up the corrections of the solver if they are at different orders.
  // The CRK hydro has already evaluated the moments, so we get to just reuse those.
  const FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CRKSPH, 0.0);
  const FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CRKSPH, Vector::zero);
  const FieldList<Dimension, SymTensor> m2 = state.fields(HydroFieldNames::m2_CRKSPH, SymTensor::zero);
  const FieldList<Dimension, ThirdRankTensor> m3 = state.fields(HydroFieldNames::m3_CRKSPH, ThirdRankTensor::zero);
  const FieldList<Dimension, FourthRankTensor> m4 = state.fields(HydroFieldNames::m4_CRKSPH, FourthRankTensor::zero);
  const FieldList<Dimension, Vector> gradm0 = state.fields(HydroFieldNames::gradM0_CRKSPH, Vector::zero);
  const FieldList<Dimension, Tensor> gradm1 = state.fields(HydroFieldNames::gradM1_CRKSPH, Tensor::zero);
  const FieldList<Dimension, ThirdRankTensor> gradm2 = state.fields(HydroFieldNames::gradM2_CRKSPH, ThirdRankTensor::zero);
  const FieldList<Dimension, FourthRankTensor> gradm3 = state.fields(HydroFieldNames::gradM3_CRKSPH, FourthRankTensor::zero);
  const FieldList<Dimension, FifthRankTensor> gradm4 = state.fields(HydroFieldNames::gradM4_CRKSPH, FifthRankTensor::zero);
  FieldList<Dimension, Scalar> QA = dataBase.newFluidFieldList(0.0, "Q A");
  FieldList<Dimension, Vector> QB;
  FieldList<Dimension, Tensor> QC;
  FieldList<Dimension, Vector> QgradA = dataBase.newFluidFieldList(Vector::zero, "Q grad A");
  FieldList<Dimension, Tensor> QgradB;
  FieldList<Dimension, ThirdRankTensor> QgradC;
  if (correctionOrder == CRKSPHSpace::LinearOrder or correctionOrder == CRKSPHSpace::QuadraticOrder) {
    QB = dataBase.newFluidFieldList(Vector::zero, "Q B");
    QgradB = dataBase.newFluidFieldList(Tensor::zero, "Q grad B");
  }
  if (correctionOrder == CRKSPHSpace::QuadraticOrder) {
    QC = dataBase.newFluidFieldList(Tensor::zero, "Q C");
    QgradC = dataBase.newFluidFieldList(ThirdRankTensor::zero, "Q grad C");
  }

  FieldList<Dimension, Tensor> QDvDx  = dataBase.newFluidFieldList(Tensor::zero, "Q Velocity Gradient");

  // Change CRKSPH weights here if need be!
  // const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  // const FieldList<Dimension, Scalar> vol = mass/massDensity;
  const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  CRKSPHSpace::computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, correctionOrder, QA, QB, QC, QgradA, QgradB, QgradC);
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr < boundaryEnd;
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(QA);
    (*boundItr)->applyFieldListGhostBoundary(QB);
    (*boundItr)->applyFieldListGhostBoundary(QC);
    (*boundItr)->applyFieldListGhostBoundary(QgradA);
    (*boundItr)->applyFieldListGhostBoundary(QgradB);
    (*boundItr)->applyFieldListGhostBoundary(QgradC);
  }
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Get the velocity gradient.
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  QDvDx.assignFields(CRKSPHSpace::gradientCRKSPH(velocity, position, vol, H, QA, QB, QC, QgradA, QgradB, QgradC, connectivityMap, ArtificialViscosity<Dimension>::QcorrectionOrder(), W, NodeCoupling()));
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr < boundaryEnd;
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(QDvDx);
  }
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Cache pointers to the velocity gradient.
  mGradVel = QDvDx;

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
  const Tensor& DvDxi = mGradVel(nodeListi, i);
  const Tensor& DvDxj = mGradVel(nodeListj, j);

  // Grab the FieldLists scaling the coefficients.
  // These incorporate things like the Balsara shearing switch or Morris & Monaghan time evolved
  // coefficients.
  const Scalar fCli = this->mClMultiplier(nodeListi, i);
  const Scalar fCqi = this->mCqMultiplier(nodeListi, i);
  const Scalar fClj = this->mClMultiplier(nodeListj, j);
  const Scalar fCqj = this->mCqMultiplier(nodeListj, j);

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
  const Vector xij = 0.5*(xi - xj);
  const Scalar gradi = (DvDxi.dot(xij)).dot(xij);
  const Scalar gradj = (DvDxj.dot(xij)).dot(xij);
  const Scalar ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  const Scalar rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  CHECK(min(ri, rj) <= 1.0);
  //const Scalar phi = limiterMM(min(ri, rj));
  Scalar phi = limiterVL(min(ri, rj));

  // If the points are getting too close, we let the Q come back full force.
  const Scalar etaij = min(etai.magnitude(), etaj.magnitude());
  // phi *= (etaij2 < etaCrit2 ? 0.0 : 1.0);
  // phi *= min(1.0, etaij2*etaij2/(etaCrit2etaCrit2));
  if (etaij < mEtaCrit) {
    phi *= exp(-FastMath::square((etaij - mEtaCrit)/mEtaFold));
  }

  // "Mike" method.
  const Vector vi1 = vi - phi*DvDxi*xij;
  const Vector vj1 = vj + phi*DvDxj*xij;
  
  const Vector vij = vi1 - vj1;
  
  // Compute mu.
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

  // The artificial internal energy.
  const Scalar ei = fshear*(-Cl*fCli*csi*(linearInExp    ? mui                : min(0.0, mui)) +
                             Cq*fCqi    *(quadInExp      ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui))));
  const Scalar ej = fshear*(-Cl*fClj*csj*(linearInExp    ? muj                : min(0.0, muj)) +
                             Cq*fCqj    *(quadInExp      ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj))));
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
CRKSPHMonaghanGingoldViscosity<Dimension>::
etaCritFrac() const {
  return mEtaCritFrac;
}

template<typename Dimension>
void
CRKSPHMonaghanGingoldViscosity<Dimension>::
etaCritFrac(const double val) {
  VERIFY(val >= 0.0);
  mEtaCritFrac = val;
}

//------------------------------------------------------------------------------
// etaFoldFrac
//------------------------------------------------------------------------------
template<typename Dimension>
double
CRKSPHMonaghanGingoldViscosity<Dimension>::
etaFoldFrac() const {
  return mEtaFoldFrac;
}

template<typename Dimension>
void
CRKSPHMonaghanGingoldViscosity<Dimension>::
etaFoldFrac(const double val) {
  VERIFY(val > 0.0);
  mEtaFoldFrac = val;
}

}
}
