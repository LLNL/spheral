//---------------------------------Spheral++----------------------------------//
// The plain old von Neuman/Richtmyer viscosity, using a MASH formalism to
// evaluate the divergence of the velocity.
//
// Created by JMO, Wed Sep 21 20:57:00 PDT 2005
//----------------------------------------------------------------------------//
#include "MASHVonNeumanViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

using std::abs;
using std::min;
using std::max;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MASHVonNeumanViscosity<Dimension>::
MASHVonNeumanViscosity(const Scalar Clinear,
                       const Scalar Cquadratic):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic),
  mVelocityDivergence(FieldList<Dimension, Scalar>::Copy),
  mCorrection(FieldList<Dimension, Tensor>::Copy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MASHVonNeumanViscosity<Dimension>::
~MASHVonNeumanViscosity() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MASHVonNeumanViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
	   const typename Dimension::Scalar time,
	   const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Prepare the FieldLists.
  mVelocityDivergence = dataBase.newFluidFieldList(Scalar(), "velocity divergence");
  mCorrection = dataBase.newFluidFieldList(Tensor(), "correction");
  FieldList<Dimension, int> flagNodeDone = dataBase.newFluidFieldList(int(), "nodes completed");
  flagNodeDone = 0;

  // Get the fluid state.
  const FieldList<Dimension, Vector> r = dataBase.fluidPosition();
  const FieldList<Dimension, Vector> v = dataBase.fluidVelocity();
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  const FieldList<Dimension, Scalar> weight = dataBase.fluidWeight();

  // First we have to determine the MASH correction tensors.
  // Iterate over the internal nodes.
  for (InternalNodeIterator<Dimension> itr = dataBase.fluidInternalNodeBegin();
       itr != dataBase.fluidInternalNodeEnd();
       ++itr) {

    // Has this node been completed yet?
    if (flagNodeDone(itr) == 0) {

      // Set the master/coarse neighbor info.
      dataBase.setMasterNeighborNodeLists(r(itr), H(itr));

      // Iterate over the masters.
      for (MasterNodeIterator<Dimension> nodeI = dataBase.fluidMasterNodeBegin();
           nodeI != dataBase.fluidMasterNodeEnd();
           ++nodeI) {
        CHECK(flagNodeDone(nodeI) == 0);

        // State for this node.
        const Vector& ri = r(nodeI);
        const SymTensor& Hi = H(nodeI);
        const Scalar wi = weight(nodeI);

        // Select the neighbor nodes.
        dataBase.setRefineNeighborNodeLists(ri, Hi);

        // Iterate over the neighbors.
        for (RefineNodeIterator<Dimension> nodeJ = r.refineNodeBegin();
             nodeJ != r.refineNodeEnd();
             ++nodeJ) {

          // Neighbor state.
          const Vector& rj = r(nodeJ);
          const SymTensor& Hj = H(nodeJ);
          const Scalar wj = weight(nodeJ);

          // Increment the normalization and correction terms.
          const Vector rij = ri - rj;

          const Vector etai = Hj*rij;
          const Vector gradWi = Hi*etai.unitVector()*W.grad(etai.magnitude(), 1.0);

          const Vector etaj = Hj*rij;
          const Vector gradWj = Hj*etaj.unitVector()*W.grad(etaj.magnitude(), 1.0);

          mCorrection(nodeI) -= (wi + wj)*rij*(gradWi + gradWj);

        }

        // Complete the gradient correction.
        CHECK(mCorrection(nodeI).Determinant() != 0.0);
        const Tensor Ci = mCorrection(nodeI).Inverse();
        mCorrection(nodeI) = Ci;
        flagNodeDone(nodeI) = 1;

      }
    }
  }

  // Now evaluate the velocity divergence.
  // Iterate over the internal nodes.
  flagNodeDone = 0;
  for (InternalNodeIterator<Dimension> itr = dataBase.fluidInternalNodeBegin();
       itr != dataBase.fluidInternalNodeEnd();
       ++itr) {

    // Has this node been completed yet?
    if (flagNodeDone(itr) == 0) {

      // Set the master/coarse neighbor info.
      dataBase.setMasterNeighborNodeLists(r(itr), H(itr));

      // Iterate over the masters.
      for (MasterNodeIterator<Dimension> nodeI = dataBase.fluidMasterNodeBegin();
           nodeI != dataBase.fluidMasterNodeEnd();
           ++nodeI) {
        CHECK(flagNodeDone(nodeI) == 0);

        // State for this node.
        const Vector& ri = r(nodeI);
        const Vector& vi = v(nodeI);
        const SymTensor& Hi = H(nodeI);
        const Scalar wi = weight(nodeI);
        const Tensor& Ci = mCorrection(nodeI);

        // Select the neighbor nodes.
        dataBase.setRefineNeighborNodeLists(ri, Hi);

        // Iterate over the neighbors.
        for (RefineNodeIterator<Dimension> nodeJ = r.refineNodeBegin();
             nodeJ != r.refineNodeEnd();
             ++nodeJ) {

          // Neighbor state.
          const Vector& rj = r(nodeJ);
          const Vector& vj = v(nodeJ);
          const SymTensor& Hj = H(nodeJ);
          const Scalar wj = weight(nodeJ);
          const Tensor& Cj = mCorrection(nodeJ);

          // Increment the velocity divergence.
          const Vector rij = ri - rj;

          const Vector etai = Hj*rij;
          const Vector gradWi = Hi*etai.unitVector()*W.grad(etai.magnitude(), 1.0);

          const Vector etaj = Hj*rij;
          const Vector gradWj = Hj*etaj.unitVector()*W.grad(etaj.magnitude(), 1.0);
          
          const Vector vji = vj - vi;
          mVelocityDivergence(nodeI) += (wi + wj)*(Ci*vji.dyad(gradWi + gradWj)).Trace();

        }

        // Flag this node as completed.
        flagNodeDone(nodeI) = 1;

      }
    }
  }

  // Finally apply the boundary conditions to the velocity divergence.
  for (ConstBoundaryIterator boundaryItr = boundaryBegin;
       boundaryItr != boundaryEnd;
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mVelocityDivergence);
  }

}

//------------------------------------------------------------------------------
// Method to calculate and return the viscous acceleration, work, and pressure,
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
template<typename Dimension>
void
MASHVonNeumanViscosity<Dimension>::
viscousEffects(typename Dimension::Vector& acceleration,
               typename Dimension::Scalar& work,
               typename Dimension::Scalar& pressure,
               const NodeIteratorBase<Dimension>& nodeI,
               const NodeIteratorBase<Dimension>& nodeJ,
               const typename Dimension::Vector& rij, 
               const typename Dimension::Vector& rijUnit, 
               const typename Dimension::Vector& vi,
               const typename Dimension::Vector& vj,
               const typename Dimension::Vector& etai,
               const typename Dimension::Vector& etaj,
               const typename Dimension::Scalar ci,
               const typename Dimension::Scalar cj,
               const typename Dimension::Scalar Pi,
               const typename Dimension::Scalar Pj,
               const typename Dimension::Scalar rhoi,
               const typename Dimension::Scalar rhoj,
               const typename Dimension::Scalar hi,
               const typename Dimension::Scalar hj,
               const typename Dimension::Vector& gradW) const {

  REQUIRE(mVelocityDivergence.fieldForNodeList(*nodeI.nodeListPtr()) != mVelocityDivergence.end());
  REQUIRE(mVelocityDivergence.fieldForNodeList(*nodeJ.nodeListPtr()) != mVelocityDivergence.end());
  REQUIRE(rhoi > 0.0);

  // Compute the viscous internal energy in terms of the velocity divergence.
  const Scalar alpha = -(this->Cl());
  const Scalar beta = this->Cq();
  const Scalar mui = min(0.0, mVelocityDivergence(nodeI));
  const Scalar Qepsi = (alpha*ci + beta*hi*mui)*hi*mui;
  CHECK(Qepsi >= 0.0);

  // Now compute the return values.
  const Vector vij = vi - vj;
  Scalar QPi = Qepsi/rhoi;
  acceleration = QPi*gradW;
  work = QPi*vij.dot(gradW);
  pressure = rhoi*Qepsi;
}

//------------------------------------------------------------------------------
// Dump the current state of the Q to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MASHVonNeumanViscosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  ArtificialViscosity<Dimension>::dumpState(file, pathName);
  file.write(mVelocityDivergence, pathName + "/velocityDivergence");
  file.write(mCorrection, pathName + "/correction");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MASHVonNeumanViscosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  ArtificialViscosity<Dimension>::restoreState(file, pathName);
  file.read(mVelocityDivergence, pathName + "/velocityDivergence");
  file.read(mCorrection, pathName + "/correction");
}  

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class MASHVonNeumanViscosity< Dim<1> >;
  template class MASHVonNeumanViscosity< Dim<2> >;
  template class MASHVonNeumanViscosity< Dim<3> >;
}
