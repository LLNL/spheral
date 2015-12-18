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
  mQuadraticInExpansion(quadraticInExpansion),
  mGradVel(FieldSpace::Reference) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldViscosity<Dimension>::
~MonaghanGingoldViscosity() {
}


//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MonaghanGingoldViscosity<Dimension>::
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
MonaghanGingoldViscosity<Dimension>::
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
  const bool balsaraShearCorrection = this->mBalsaraShearCorrection;
  const FieldSpace::FieldList<Dimension, Scalar>& rvAlphaQ = this->reducingViscosityMultiplierQ();
  const FieldSpace::FieldList<Dimension, Scalar>& rvAlphaL = this->reducingViscosityMultiplierL();

  // Are we applying the shear corrections?
  //const Scalar fsheari = (balsaraShearCorrection ? this->mShearMultiplier(nodeListi, i) : 1.0);
  //const Scalar fshearj = (balsaraShearCorrection ? this->mShearMultiplier(nodeListj, j) : 1.0);

  Scalar fshear = 1.0;
  Scalar fsheari = fshear;
  Scalar fshearj = fshear;
  const Tensor& DvDxi = mGradVel(nodeListi, i);
  const Tensor& DvDxj = mGradVel(nodeListj, j);
  if (balsaraShearCorrection) {
    const Scalar csneg = this->negligibleSoundSpeed();
    const Scalar hiinv = Hi.Trace()/Dimension::nDim;
    const Scalar hjinv = Hj.Trace()/Dimension::nDim;
    const Scalar ci = max(csneg, csi);
    const Scalar cj = max(csneg, csj);
    const Scalar fi = abs(DvDxi.Trace())/(this->curlVelocityMagnitude(DvDxi) + abs(DvDxi.Trace()) + eps2*ci*hiinv);
    const Scalar fj = abs(DvDxj.Trace())/(this->curlVelocityMagnitude(DvDxj) + abs(DvDxj.Trace()) + eps2*cj*hjinv);
    fshear = min(fi, fj);
    //fsheari = fi;
    //fshearj = fj;
    fsheari = fshear;
    fshearj = fshear;
  }

  // Compute mu.
  const Vector vij = vi - vj;
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

  // The artificial internal energy.
  // const Scalar ei = fshear*(-Cl*csi*(mLinearInExpansion    ? mui                : min(0.0, mui)) +
  //                            Cq     *(mQuadraticInExpansion ? -sgn(mui)*mui*mui : FastMath::square(min(0.0, mui))));
  // const Scalar ej = fshear*(-Cl*csj*(mLinearInExpansion    ? muj                : min(0.0, muj)) +
  //                            Cq     *(mQuadraticInExpansion ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj))));
  const Scalar ei = fsheari*(-Cl*rvAlphaL(nodeListi,i)*csi*(mLinearInExpansion    ? mui                : min(0.0, mui)) +
                              Cq *rvAlphaQ(nodeListi,i)   *(mQuadraticInExpansion ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui)))) ;
  const Scalar ej = fshearj*(-Cl*rvAlphaL(nodeListj,j)*csj*(mLinearInExpansion    ? muj                : min(0.0, muj)) +
                              Cq *rvAlphaQ(nodeListj,j)    *(mQuadraticInExpansion ? -sgn(muj)*muj*muj : FastMath::square(min(0.0, muj))));
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
}
