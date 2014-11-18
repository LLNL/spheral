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
  sigma = CSPHSpace::gradientCSPH(velocity, position, vol, H, A, B, C, D, gradA, gradB, connectivityMap, W);

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

  // Apply boundary conditions to div velocity.
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
  }
  // for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
  //      boundItr != boundaryEnd;
  //      ++boundItr) {
  //   (*boundItr)->finalizeGhostBoundary();
  // }
}

}
}
