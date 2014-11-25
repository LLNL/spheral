//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- The base class for all ArtificialViscosities in 
// Spheral++.
//----------------------------------------------------------------------------//

#include <algorithm>

#include "TensorCSPHViscosity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"
#include "FileIO/FileIO.hh"
#include "CSPH/gradientCSPH.hh"

using namespace std;

namespace Spheral {
namespace ArtificialViscositySpace {

using std::vector;
using NodeSpace::NodeList;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using BoundarySpace::Boundary;
using KernelSpace::TableKernel;
using NeighborSpace::Neighbor;
using FileIOSpace::FileIO;
using NeighborSpace::ConnectivityMap;
using NodeSpace::FluidNodeList;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorCSPHViscosity<Dimension>::
TensorCSPHViscosity(Scalar Clinear, Scalar Cquadratic):
  mGradVel(FieldSpace::Copy),
  TensorMonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorCSPHViscosity<Dimension>::
~TensorCSPHViscosity() {
}

//------------------------------------------------------------------------------
// Method to apply the viscous acceleration, work, and pressure, to the derivatives
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
     typename Dimension::Tensor>
TensorCSPHViscosity<Dimension>::
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
    const double Cl = this->mClinear;
    const double Cq = this->mCquadratic;
    const double eps2 = this->mEpsilon2;
    const bool balsara = this->mBalsaraShearCorrection;
    const bool limiter = this->mLimiterSwitch;

    // State of node I.
    const Tensor& _sigmai = this->mSigma(nodeListi, i);
    const Scalar fsheari = (balsara ?
                            this->mShearMultiplier(nodeListi, i) :
                            1.0);
    CHECK(fsheari >= 0.0 && fsheari <= 1.0);

    // State for node J.
    const Tensor& _sigmaj = this->mSigma(nodeListj, j);
    const Scalar fshearj = (balsara ?
                            this->mShearMultiplier(nodeListj, j) :
                            1.0);
    CHECK(fshearj >= 0.0 && fshearj <= 1.0);

    // Some more geometry.
    const Scalar xij2 = xij.magnitude2();
    const Vector xijUnit = xij.unitVector();
    const Scalar hi2 = xij2/(etai.magnitude2() + tiny);
    const Scalar hj2 = xij2/(etaj.magnitude2() + tiny);
    const Scalar hi = sqrt(hi2);
    const Scalar hj = sqrt(hj2);

    // BOOGA!
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
    Tensor Qepsi = fsheari*(-Cl*csi*mui.Transpose() + Cq*mui*mui);
    if (limiter) Qepsi = this->calculateLimiter(vi, vj, csi, csj, hi, hj, nodeListi, i)*Qepsi;

    const Tensor muj = hj*sigmaj;
    Tensor Qepsj = fshearj*(-Cl*csj*muj.Transpose() + Cq*muj*muj);
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
TensorCSPHViscosity<Dimension>::
calculateSigmaAndGradDivV(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const TableKernel<Dimension>& W,
                          typename TensorCSPHViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                          typename TensorCSPHViscosity<Dimension>::ConstBoundaryIterator boundaryEnd) {

  FieldList<Dimension, Tensor>& sigma = ArtificialViscosity<Dimension>::mSigma;
  FieldList<Dimension, Vector>& gradDivVelocity = ArtificialViscosity<Dimension>::mGradDivVelocity;

  // Get the necessary state.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  const FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  const FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);

  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const int numNodeLists = dataBase.numFluidNodeLists();

  // Compute the basic velocity gradient.
  const FieldList<Dimension, Scalar> vol = mass/rho;
  mGradVel = CSPHSpace::gradientCSPH(velocity, position, vol, H, A, B, C, D, gradA, gradB, connectivityMap, W);
  sigma = mGradVel;
  sigma.copyFields();

  // Compute sigma and build the velocity divergence.
  FieldList<Dimension, Scalar> divVel = dataBase.newFluidFieldList(0.0, "velocity divergence");
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;
      Tensor& sigmai = sigma(nodeListi, i);

      // Update the velocity divergence.
      divVel(nodeListi, i) = sigmai.Trace();

      // Now limit to just negative eigen-values.  This is 'cause we only
      // care about convergent geometries for the Q.
      const SymTensor sigmai_s = sigmai.Symmetric();
      const Tensor sigmai_a = sigmai.SkewSymmetric();
      typename SymTensor::EigenStructType eigeni = sigmai_s.eigenVectors();
      sigmai = constructTensorWithMinDiagonal(eigeni.eigenValues, 0.0);
      sigmai.rotationalTransform(eigeni.eigenVectors);
      // sigmai += sigmai_a;
    }
  }

  // Apply boundary conditions.
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr < boundaryEnd;
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(divVel);
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Compute the gradient of div vel.
  gradDivVelocity = CSPHSpace::gradientCSPH(divVel, position, vol, H, A, B, C, D, gradA, gradB, connectivityMap, W);

  // Apply boundary conditions.
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr < boundaryEnd;
       ++boundItr) {
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
}
