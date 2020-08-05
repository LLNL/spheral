//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor 
// formalism.
//----------------------------------------------------------------------------//
#include "TensorMonaghanGingoldViscosity.hh"
#include "DataOutput/Restart.hh"
#include "Boundary/Boundary.hh"
#include "Geometry/EigenStruct.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"

#include "Utilities/DBC.hh"

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
TensorMonaghanGingoldViscosity<Dimension>::
TensorMonaghanGingoldViscosity(Scalar Clinear, Scalar Cquadratic):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic) {
  this->calculateSigma(true);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorMonaghanGingoldViscosity<Dimension>::
~TensorMonaghanGingoldViscosity() {
}

//------------------------------------------------------------------------------
// Method to apply the viscous acceleration, work, and pressure, to the derivatives
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
     typename Dimension::Tensor>
TensorMonaghanGingoldViscosity<Dimension>::
Piij(const unsigned nodeListi, const unsigned i, 
     const unsigned nodeListj, const unsigned j,
     const Vector& xi,
     const Vector& etai,
     const Vector& vi,
     const Scalar rhoi,
     const Scalar csi,
     const SymTensor& /*Hi*/,
     const Vector& xj,
     const Vector& etaj,
     const Vector& vj,
     const Scalar rhoj,
     const Scalar csj,
     const SymTensor& /*Hj*/) const {

  // If the nodes are not closing, then skip the rest and the Q is zero.
  const Vector xij = xi - xj;
  const Vector vij = vi - vj;
  if (vij.dot(xij) < 0.0) {
          
    const double tiny = 1.0e-20;
    double Cl = this->mClinear;
    double Cq = this->mCquadratic;
    const double eps2 = this->mEpsilon2;
    const bool limiter = this->mLimiterSwitch;

    // Grab the FieldLists scaling the coefficients.
    // These incorporate things like the Balsara shearing switch or Morris & Monaghan time evolved
    // coefficients.
    const Scalar fCli = this->mClMultiplier(nodeListi, i);
    const Scalar fCqi = this->mCqMultiplier(nodeListi, i);
    const Scalar fClj = this->mClMultiplier(nodeListj, j);
    const Scalar fCqj = this->mCqMultiplier(nodeListj, j);
    const Scalar fshear = std::max(this->mShearCorrection(nodeListi, i), this->mShearCorrection(nodeListj, j));
    Cl *= 0.5*(fCli + fClj)*fshear;
    Cq *= 0.5*(fCqi + fCqj)*fshear;

    // Some more geometry.
    const Scalar xij2 = xij.magnitude2();
    const Vector xijUnit = xij.unitVector();
    const Scalar hi2 = xij2/(etai.magnitude2() + tiny);
    const Scalar hj2 = xij2/(etaj.magnitude2() + tiny);
    const Scalar hi = sqrt(hi2);
    const Scalar hj = sqrt(hj2);

    // BOOGA!
    const Tensor& _sigmai = this->mSigma(nodeListi, i);
    const Tensor& _sigmaj = this->mSigma(nodeListj, j);
    Tensor sigmai = _sigmai;
    Tensor sigmaj = _sigmaj;
    {
      const Tensor R = rotationMatrix(xijUnit);
      const Tensor Rinverse = R.Transpose();
      const Vector thpt1 = sqrt(xij2)*(R*vij);
      const Vector deltaSigmai = thpt1/(xij2 + eps2*hi2);
      const Vector deltaSigmaj = thpt1/(xij2 + eps2*hj2);
      sigmai.rotationalTransform(R);
      sigmaj.rotationalTransform(R);
      sigmai.setColumn(0, deltaSigmai);
      sigmaj.setColumn(0, deltaSigmaj);
      sigmai.rotationalTransform(Rinverse);
      sigmaj.rotationalTransform(Rinverse);
    }
    // BOOGA!

    // Calculate the tensor viscous internal energy.
    const Tensor mui = hi*sigmai;
    Tensor Qepsi = -Cl*csi*mui.Transpose() + Cq*mui*mui;
    if (limiter) Qepsi = this->calculateLimiter(vi, vj, csi, csj, hi, hj, nodeListi, i)*Qepsi;

    const Tensor muj = hj*sigmaj;
    Tensor Qepsj = -Cl*csj*muj.Transpose() + Cq*muj*muj;
    if (limiter) Qepsj = this->calculateLimiter(vj, vi, csj, csi, hj, hi, nodeListj, j)*Qepsj;

    // We now have enough to compute Pi!
    const Tensor QPii = Qepsi/rhoi;
    const Tensor QPij = Qepsj/rhoj;
    return make_pair(QPii, QPij);

  } else {

    return make_pair(Tensor::zero, Tensor::zero);

  }
}

}

