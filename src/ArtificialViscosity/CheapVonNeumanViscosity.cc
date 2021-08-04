//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "CheapVonNeumanViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CheapVonNeumanViscosity<Dimension>::
CheapVonNeumanViscosity():
  VonNeumanViscosity<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
CheapVonNeumanViscosity<Dimension>::
CheapVonNeumanViscosity(Scalar Clinear, Scalar Cquadratic):
  VonNeumanViscosity<Dimension>(Clinear, Cquadratic) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CheapVonNeumanViscosity<Dimension>::
~CheapVonNeumanViscosity() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CheapVonNeumanViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
           const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Verify that the internal pressure field is properly initialized for this
  // set of fluid node lists.
  if (mViscousEnergy.numFields() != dataBase.numFluidNodeLists()) {
    FieldList<Dimension, Scalar> thpt(FieldList<Dimension, Scalar>::Copy);
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator fluidNodeListItr =
           dataBase.fluidNodeListBegin();
         fluidNodeListItr < dataBase.fluidNodeListEnd();
         ++fluidNodeListItr) {
      thpt.appendField(Field<Dimension, Scalar>(**fluidNodeListItr));
    }
    mViscousEnergy = thpt;
  }
  CHECK(mViscousEnergy.numFields() == dataBase.numFluidNodeLists());

  // Get the fluid state.
  const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  const FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  const FieldList<Dimension, Scalar> specificEnergy = dataBase.fluidSpecificThermalEnergy();
  const FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  const FieldList<Dimension, Scalar> weight = dataBase.fluidWeight();
  const FieldList<Dimension, Scalar> pressure = dataBase.fluidPressure();
  const FieldList<Dimension, Scalar> soundSpeed = dataBase.fluidSoundSpeed();
  const FieldList<Dimension, Tensor> DvDx = dataBase.fluidDvelocityDx();

  // Set the viscous energy for each fluid node.
  const Scalar alpha = -Cl();
  const Scalar beta = Cq();
  for (NodeIDIterator<Dimension> nodeItr = dataBase.neighborInternalNodeBegin();
       nodeItr < dataBase.neighborInternalNodeEnd();
       ++nodeItr) {
    const Scalar mui = min(0.0, DvDx(nodeItr).Trace());
    const Scalar l = nodeItr.nodeListPtr()->neighborPtr()->nodeExtentField()(nodeItr).maxElement();
    mViscousEnergy(nodeItr) = (alpha*soundSpeed(nodeItr) + beta*l*mui)*l*mui;
  }

  // Finally apply the boundary conditions to the viscous energy FieldList.
  for (ConstBoundaryIterator boundaryItr = boundaryBegin;
       boundaryItr < boundaryEnd;
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mViscousEnergy);
  }

}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class CheapVonNeumanViscosity< Dim<1> >;
  template class CheapVonNeumanViscosity< Dim<2> >;
  template class CheapVonNeumanViscosity< Dim<3> >;
}
