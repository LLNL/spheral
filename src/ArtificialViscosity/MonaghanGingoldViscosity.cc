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
  FieldSpace::FieldList<Dimension, Scalar>& rvAlpha = this->reducingViscosityMultiplier();

  // Are we applying the shear corrections?
  const Vector vij = vi - vj;
  Scalar fshear = 1.0;

    //add scalar for rvAlpha grab and use later
    
  if (balsaraShearCorrection) {
    fshear = abs(vij.unitVector().dot(vij.unitVector()));

//     const Scalar hmax = max(Hi.Inverse().Trace(), Hj.Inverse().Trace())/Dimension::nDim;
//     const Scalar hmaxinv = 1.0/hmax; 
//     const Vector xij = xi - xj;
//     Tensor DvDx;
//     for (unsigned i = 0; i != Dimension::nDim; ++i) {
//       for (unsigned j = 0; j != Dimension::nDim; ++j) {
//         DvDx(i,j) = vij(i)/max(xij(j), 0.1*hmax);
//       }
//     }
//     const Scalar div = abs(DvDx.Trace());
//     const Scalar curl = curlVelocityMagnitude(DvDx);
//     const Scalar cs = max(this->negligibleSoundSpeed(), max(csi, csj));
//     fshear = div/(div + curl + this->epsilon2()*cs*hmaxinv);
  }

//   const Scalar fsheari = (this->mBalsaraShearCorrection ? 
//                           this->mShearMultiplier(nodeListi, i) :
//                           1.0);
//   const Scalar fshearj = (this->mBalsaraShearCorrection ? 
//                           this->mShearMultiplier(nodeListj, j) :
//                           1.0);

  // Compute mu.
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

  // The artificial internal energy.
  const Scalar ei = fshear*(-Cl*csi*(mLinearInExpansion    ? mui                : min(0.0, mui)) +
                            Cq     *(mQuadraticInExpansion ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui))))*rvAlpha(nodeListi,i);
  const Scalar ej = fshear*(-Cl*csj*(mLinearInExpansion    ? muj                : min(0.0, muj)) +
                            Cq     *(mQuadraticInExpansion ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj))))*rvAlpha(nodeListj,j);
  CHECK(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion));
  CHECK(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion));

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
