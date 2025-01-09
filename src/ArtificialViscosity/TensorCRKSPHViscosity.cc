//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- The base class for all ArtificialViscosities in 
// Spheral++.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"
#include "RK/RKFieldNames.hh"
#include "RK/gradientRK.hh"

#include "TensorCRKSPHViscosity.hh"

#include <algorithm>

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
TensorCRKSPHViscosity<Dimension>::
TensorCRKSPHViscosity(Scalar Clinear, Scalar Cquadratic):
  TensorMonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic),
  mGradVel(FieldStorageType::CopyFields){
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorCRKSPHViscosity<Dimension>::
~TensorCRKSPHViscosity() {
}

//------------------------------------------------------------------------------
// Method to apply the viscous acceleration, work, and pressure, to the derivatives
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
     typename Dimension::Tensor>
TensorCRKSPHViscosity<Dimension>::
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

  // Estimate the velocity difference.
  Vector vij = vi - vj;
  const Vector xij = xi - xj;
  const Tensor& DvDxi = mGradVel(nodeListi, i);
  const Tensor& DvDxj = mGradVel(nodeListj, j);
  const Vector vi1 = vi - 0.5*(DvDxi.dot(xij));
  const Vector vj1 = vj + 0.5*(DvDxj.dot(xij));
  const Vector vij1 = vi1 - vj1;
  if (((vi1 - vj).dot(vij) > 0.0) and
      ((vi - vj1).dot(vij) > 0.0) and
      (vij1.dot(vij) > 0.0)) {
    const Vector vijhat = vij.unitVector();
    // const bool barf = (vij.magnitude() > 1.0e-5);
    // if (barf) cerr << "Cutting vij from " << vij << " to ";
    vij = min(vij.magnitude(), vij1.magnitude())*vijhat;
    // if (barf) cerr << vij << endl;
  }

  // If the nodes are not closing, then skip the rest and the Q is zero.
  if (vij.dot(xij) < 0.0) {
          
    const double tiny = 1.0e-20;
    double Cl = this->mClinear;
    double Cq = this->mCquadratic;
    const double eps2 = this->mEpsilon2;
    //const bool balsara = this->mBalsaraShearCorrection;
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

//------------------------------------------------------------------------------
// Compute the internal background sigma and grad-div-v fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorCRKSPHViscosity<Dimension>::
calculateSigmaAndGradDivV(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& /*derivs*/,
                          const TableKernel<Dimension>& /*W*/,
                          typename TensorCRKSPHViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                          typename TensorCRKSPHViscosity<Dimension>::ConstBoundaryIterator boundaryEnd) {

  const auto order = ArtificialViscosity<Dimension>::QcorrectionOrder();
  auto& sigma = ArtificialViscosity<Dimension>::mSigma;
  auto& gradDivVelocity = ArtificialViscosity<Dimension>::mGradDivVelocity;

  // Get the necessary state.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto WR = state.template get<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(order));
  const auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());

  const auto& connectivityMap = dataBase.connectivityMap();
  const auto numNodeLists = dataBase.numFluidNodeLists();

  // Compute the basic velocity gradient.
  const auto vol = mass/rho;
  mGradVel = gradientRK(velocity, position, vol, H, connectivityMap, WR, corrections,  NodeCoupling());
  sigma = mGradVel;
  sigma.copyFields();

  // Compute sigma and build the velocity divergence.
  auto divVel = dataBase.newFluidFieldList(0.0, "velocity divergence");
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
    for (auto iItr = connectivityMap.begin(nodeListi); iItr < connectivityMap.end(nodeListi); ++iItr) {
      const auto i = *iItr;
      auto& sigmai = sigma(nodeListi, i);

      // Update the velocity divergence.
      divVel(nodeListi, i) = sigmai.Trace();

      // Now limit to just negative eigen-values.  This is 'cause we only
      // care about convergent geometries for the Q.
      const auto sigmai_s = sigmai.Symmetric();
      const auto sigmai_a = sigmai.SkewSymmetric();
      auto eigeni = sigmai_s.eigenVectors();
      sigmai = constructTensorWithMinDiagonal(eigeni.eigenValues, 0.0);
      sigmai.rotationalTransform(eigeni.eigenVectors);
      // sigmai += sigmai_a;
    }
  }

  // Apply boundary conditions.
  for (auto boundItr = boundaryBegin; boundItr < boundaryEnd; ++boundItr) (*boundItr)->applyFieldListGhostBoundary(divVel);
  for (auto boundItr = boundaryBegin; boundItr < boundaryEnd; ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Compute the gradient of div vel.
  gradDivVelocity = gradientRK(divVel, position, vol, H, connectivityMap, WR, corrections, NodeCoupling());

  // Apply boundary conditions.
  for (auto boundItr = boundaryBegin; boundItr < boundaryEnd; ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(sigma);
    (*boundItr)->applyFieldListGhostBoundary(gradDivVelocity);
    (*boundItr)->applyFieldListGhostBoundary(mGradVel);
  }
  // for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
  //      boundItr != boundaryEnd;
  //      ++boundItr) {
  //   (*boundItr)->finalizeGhostBoundary();
  // }
}

}
