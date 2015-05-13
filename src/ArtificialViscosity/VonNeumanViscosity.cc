//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "VonNeumanViscosity.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPH/gradientCRKSPH.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using KernelSpace::TableKernel;
using BoundarySpace::Boundary;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
VonNeumanViscosity<Dimension>::
VonNeumanViscosity(Scalar Clinear, Scalar Cquadratic):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic),
  mViscousEnergy(FieldSpace::Copy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VonNeumanViscosity<Dimension>::
~VonNeumanViscosity() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VonNeumanViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           const State<Dimension>& state,
           const StateDerivatives<Dimension>& derivs,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
	   const typename Dimension::Scalar time,
	   const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Make sure the CRKSPH corrections have had boundaries completed.
  for (ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Verify that the internal pressure field is properly initialized for this
  // set of fluid node lists.
  dataBase.resizeFluidFieldList(mViscousEnergy, 0.0, "viscous pressure", true);
  CHECK(mViscousEnergy.numFields() == dataBase.numFluidNodeLists());

  // Get the fluid state.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  const FieldList<Dimension, Scalar> vol = mass/massDensity;

  // Get the fluid velocity gradient.
  const FieldList<Dimension, Tensor> velocityGradient = CRKSPHSpace::gradientCRKSPH(velocity,
                                                                                position,
                                                                                vol,
                                                                                H,
                                                                                A,
                                                                                B,
                                                                                gradA,
                                                                                gradB,
                                                                                connectivityMap,
                                                                                W);

  // Set the viscous energy for each fluid node.
  const Scalar Cl = this->Cl();
  const Scalar Cq = this->Cq();
  const unsigned numNodeLists = dataBase.numFluidNodeLists();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = velocityGradient[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar mui = min(0.0, velocityGradient(nodeListi, i).Trace());
      const Scalar l = 1.0/(H(nodeListi, i).eigenValues().maxElement());
      mViscousEnergy(nodeListi, i) = (-Cl*soundSpeed(nodeListi, i) + Cq*l*mui)*l*mui;
    }
  }

  // Finally apply the boundary conditions to the viscous energy FieldList.
  for (ConstBoundaryIterator boundaryItr = boundaryBegin;
       boundaryItr < boundaryEnd;
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mViscousEnergy);
  }
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<typename Dimension>
std::pair<typename Dimension::Tensor, typename Dimension::Tensor>
VonNeumanViscosity<Dimension>::
Piij(const unsigned nodeListi, const unsigned i, 
     const unsigned nodeListj, const unsigned j,
     const typename Dimension::Vector& xi,
     const typename Dimension::Vector& etai,
     const typename Dimension::Vector& vi,
     const typename Dimension::Scalar rhoi,
     const typename Dimension::Scalar csi,
     const typename Dimension::SymTensor& Hi,
     const typename Dimension::Vector& xj,
     const typename Dimension::Vector& etaj,
     const typename Dimension::Vector& vj,
     const typename Dimension::Scalar rhoj,
     const typename Dimension::Scalar csj,
     const typename Dimension::SymTensor& Hj) const {

  REQUIRE(rhoi > 0.0);

  return make_pair(mViscousEnergy(nodeListi, i)/rhoi*Tensor::one,
                   mViscousEnergy(nodeListj, j)/rhoj*Tensor::one);
}

//------------------------------------------------------------------------------
// Access the viscous energy.
//------------------------------------------------------------------------------
template<typename Dimension>
const FieldList<Dimension, typename Dimension::Scalar>&
VonNeumanViscosity<Dimension>::
viscousEnergy() const {
  return mViscousEnergy;
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
  file.write(mViscousEnergy, pathName + "/viscousEnergy");
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
