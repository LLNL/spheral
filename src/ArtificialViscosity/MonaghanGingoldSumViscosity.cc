//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "MonaghanGingoldSumViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldSumViscosity<Dimension>::
MonaghanGingoldSumViscosity():
  MonaghanGingoldViscosity<Dimension>(),
  mViscousEnergy(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldSumViscosity<Dimension>::
MonaghanGingoldSumViscosity(Scalar Clinear, Scalar Cquadratic):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic),
  mViscousEnergy(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldSumViscosity<Dimension>::
~MonaghanGingoldSumViscosity() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MonaghanGingoldSumViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           typename vector<Boundary<Dimension>*>::const_iterator boundaryBegin,
           typename vector<Boundary<Dimension>*>::const_iterator boundaryEnd,
           const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef NodeIDIterator<Dimension> IDIterator;

  // Verify that the internal pressure field is properly initialized for this
  // set of fluid node lists.
  if (mViscousEnergy.numFields() != dataBase.numFluidNodeLists()) {
    FieldList<Dimension, Scalar> thpt(FieldStorageType::CopyFields);
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator fluidNodeListItr = 
           dataBase.fluidNodeListBegin();
         fluidNodeListItr < dataBase.fluidNodeListEnd();
         ++fluidNodeListItr) {
      thpt.appendField(Field<Dimension, Scalar>(**fluidNodeListItr));
    }
    mViscousEnergy = thpt;
  } else {
    mViscousEnergy.Zero();
  }
  CHECK(mViscousEnergy.numFields() == dataBase.numFluidNodeLists());

  // Create a list of flags per node to check when a given node is done.
  vector< vector<bool> > flagNodeDone(dataBase.numFluidNodeLists());
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator nodeListItr =
         dataBase.fluidNodeListBegin();
       nodeListItr < dataBase.fluidNodeListEnd();
       ++nodeListItr) {
    flagNodeDone[distance(dataBase.fluidNodeListBegin(), nodeListItr)].resize((*nodeListItr)->numInternalNodes(), false);
  }

  // Get the fluid state.
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  const FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  const FieldList<Dimension, Scalar> weight = dataBase.fluidWeight();
  const FieldList<Dimension, Scalar> soundSpeed = dataBase.fluidSoundSpeed();

  // Get the linear and quadratic coefficients.
  const Scalar alpha = -Cl();
  const Scalar beta = Cq();

  // Loop over the internal nodes in the DataBase.
  for (IDIterator nodeItr = dataBase.fluidInternalNodeBegin();
       nodeItr != dataBase.fluidInternalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      dataBase.setMasterNeighborNodeLists(position(nodeItr), Hfield(nodeItr));

      // Now loop over all the master nodes.
      for (IDIterator masterItr = dataBase.masterNodeBegin();
           masterItr < dataBase.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

        // Master fluid state.
        const Vector& ri = position(masterItr);
        const Vector& vi = velocity(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar weighti = weight(masterItr);
        const Scalar ci = soundSpeed(masterItr);

        // Set the refined neighbor information for this master node.
        dataBase.setRefineNeighborNodeLists(position(masterItr), Hfield(masterItr));

        // Loop over the refined neighbors for this master node and calculate
        // an estimate of del^2 v.
        Scalar mui = 0.0;
        for (IDIterator neighborItr = dataBase.refineNodeBegin();
             neighborItr < dataBase.refineNodeEnd();
             ++neighborItr) {

          // Neighbor fluid state.
          const Vector& rj = position(neighborItr);
          const Vector& vj = velocity(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar weightj = weight(neighborItr);

          // Get the kernel weighting.
          Vector rij = ri - rj;
          Vector etai = Hi*rij;
          Vector etaj = Hj*rij;
          Scalar Wij;
          switch(neighborItr.nodeListPtr()->neighborPtr()->neighborSearchType()) {

          case Neighbor<Dimension>::GatherScatter:
            Wij = 0.5*(weighti*W(Hi*rij, Hi) + weightj*W(Hj*rij, Hj));
            break;

          case Neighbor<Dimension>::Gather:
            Wij = weighti*W(Hi*rij, Hi);
            break;

          case Neighbor<Dimension>::Scatter:
            Wij = weightj*W(Hj*rij, Hj);
            break;
          }

          // This pairs contribution to the estimate of del^2 v.
          Vector vij = vi - vj;
          mui += vij.dot(rij)/(rij.magnitude2() + epsilon2())*Wij;
        }

        // Use the estimate of del^2 v to find the viscous energy for this master
        // node.
        const Scalar l = masterItr.nodeListPtr()->neighborPtr()->nodeExtentField()(masterItr).maxElement();
        mui = min(0.0, mui);
        mViscousEnergy(masterItr) = (alpha*ci+ beta*l*mui)*l*mui;

        // This master node is done.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // Finally apply the boundary conditions to the viscous energy FieldList.
  for (ConstBoundaryIterator boundaryItr = boundaryBegin;
       boundaryItr != boundaryEnd;
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mViscousEnergy);
  }
  for (ConstBoundaryIterator boundaryItr = boundaryBegin;
       boundaryItr != boundaryEnd;
       ++boundaryItr) {
    (*boundaryItr)->finalizeGhostBoundary();
  }

}

//------------------------------------------------------------------------------
// Test if the ArtificialViscosity is valid, i.e., ready to use.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MonaghanGingoldSumViscosity<Dimension>::valid() const {
  return (MonaghanGingoldViscosity<Dimension>::valid());
//           mViscousEnergy.numFields() > 0);
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class MonaghanGingoldSumViscosity< Dim<1> >;
  template class MonaghanGingoldSumViscosity< Dim<2> >;
  template class MonaghanGingoldSumViscosity< Dim<3> >;
}
