//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "VonNeumanViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListFunctions.hh"
#include "Field/FieldListFunctionsMash.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

using Spheral::DataBaseSpace::DataBase;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::KernelSpace::TableKernel;
using Spheral::BoundarySpace::Boundary;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VonNeumanViscosity<Dimension>::
VonNeumanViscosity():
  ArtificialViscosity<Dimension>(),
  mViscousEnergy(FieldList<Dimension, Scalar>::Copy) {
#ifdef DEBUG
  cerr << "VonNeumanViscosity::VonNeumanViscosity()" << endl;
#endif
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
VonNeumanViscosity<Dimension>::
VonNeumanViscosity(Scalar Clinear, Scalar Cquadratic):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic),
  mViscousEnergy(FieldList<Dimension, Scalar>::Copy) {
#ifdef DEBUG
  cerr << "VonNeumanViscosity::VonNeumanViscosity(Cl, Cq)" << endl;
#endif
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VonNeumanViscosity<Dimension>::
~VonNeumanViscosity() {
#ifdef DEBUG
  cerr << "VonNeumanViscosity::~VonNeumanViscosity()" << endl;
#endif
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VonNeumanViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
	   const typename Dimension::Scalar time,
	   const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {
#ifdef DEBUG
  cerr << "VonNeumanViscosity::initialize()" << endl;
#endif

  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Verify that the internal pressure field is properly initialized for this
  // set of fluid node lists.
  if (mViscousEnergy.numFields() != dataBase.numFluidNodeLists()) {
    FieldList<Dimension, Scalar> thpt(FieldList<Dimension, Scalar>::Copy);
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator fluidNodeListItr = dataBase.fluidNodeListBegin();
         fluidNodeListItr < dataBase.fluidNodeListEnd();
         ++fluidNodeListItr) {
      thpt.appendField(Field<Dimension, Scalar>(**fluidNodeListItr));
    }
    mViscousEnergy = thpt;
  } else {
    mViscousEnergy.Zero();
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

  // Get the fluid velocity divergence.
  const FieldList<Dimension, Scalar> velocityDivergence = 
    Spheral::FieldSpace::divergence<Dimension, Vector>(velocity,
                                                       position,
                                                       weight,
                                                       mass,
                                                       massDensity,
                                                       Hfield,
                                                       W);

  // Set the viscous energy for each fluid node.
  const Scalar alpha = -Cl();
  const Scalar beta = Cq();
  for (NodeIDIterator<Dimension> nodeItr = dataBase.internalNodeBegin();
       nodeItr < dataBase.internalNodeEnd();
       ++nodeItr) {
    const Scalar mui = min(0.0, velocityDivergence(nodeItr));
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

//------------------------------------------------------------------------------
// Method to calculate and return the viscous acceleration, work, and pressure,
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
template<typename Dimension>
void
VonNeumanViscosity<Dimension>::
viscousEffects(typename Dimension::Vector& acceleration,
               typename Dimension::Scalar& work,
               typename Dimension::Scalar& pressure,
               const NodeIDIterator<Dimension>& nodeI,
               const NodeIDIterator<Dimension>& nodeJ,
               const typename Dimension::Vector& rij, 
               const typename Dimension::Vector& vi,
               const typename Dimension::Vector& vj,
               const typename Dimension::Vector& etai,
               const typename Dimension::Vector& etaj,
               const typename Dimension::Scalar ci,
               const typename Dimension::Scalar cj,
               const typename Dimension::Scalar rhoi,
               const typename Dimension::Scalar rhoj,
               const typename Dimension::Vector& gradW) const {

  REQUIRE(rhoi > 0.0);

  // viscousInternalEnergy does most of the work.
  Vector vij = vi - vj;
  Scalar Qepsi = mViscousEnergy(nodeI);
  CHECK(Qepsi >= 0.0);
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
VonNeumanViscosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Call the ancestor method.
  ArtificialViscosity<Dimension>::dumpState(file, pathName);

  // Dump the viscous energy.
  file.write(viscousInternalEnergy(), pathName + "/viscousEnergy");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VonNeumanViscosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // Call the ancestor method.
  ArtificialViscosity<Dimension>::restoreState(file, pathName);

  // Dump the viscous energy.
  file.read(mViscousEnergy, pathName + "/viscousEnergy");
}  
}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
namespace ArtificialViscositySpace {
template class VonNeumanViscosity< Dim<1> >;
template class VonNeumanViscosity< Dim<2> >;
template class VonNeumanViscosity< Dim<3> >;
}
}
