//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "RK/RKFieldNames.hh"
#include "RK/gradientRK.hh"
#include "Hydro/HydroFieldNames.hh"

#include "VonNeumanViscosity.hh"

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
VonNeumanViscosity<Dimension>::
VonNeumanViscosity(Scalar Clinear, Scalar Cquadratic):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic),
  mViscousEnergy(FieldStorageType::CopyFields) {
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

  // Let the base class do it's thing.
  ArtificialViscosity<Dimension>::initialize(dataBase, state, derivs, boundaryBegin, boundaryEnd, time, dt, W);
  const auto order = this->QcorrectionOrder();

  // Make sure the RK corrections have had boundaries completed.
  for (auto boundItr = boundaryBegin; boundItr < boundaryEnd; ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Verify that the internal pressure field is properly initialized for this
  // set of fluid node lists.
  dataBase.resizeFluidFieldList(mViscousEnergy, 0.0, "viscous pressure", true);
  CHECK(mViscousEnergy.numFields() == dataBase.numFluidNodeLists());

  // Get the fluid state.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  specificEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto  soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto  vol = mass/massDensity;
  const auto  WR = state.template get<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(order));
  const auto  corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());

  // We'll compute the higher-accuracy RK gradient.
  const auto velocityGradient = gradientRK(velocity,
                                           position,
                                           vol,
                                           H,
                                           connectivityMap,
                                           WR,
                                           corrections,
                                           NodeCoupling());

  // Set the viscous energy for each fluid node.
  const auto Cl = this->Cl();
  const auto Cq = this->Cq();
  const auto numNodeLists = dataBase.numFluidNodeLists();
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
    const auto n = velocityGradient[nodeListi]->numInternalElements();
    for (auto i = 0u; i != n; ++i) {
      const auto mui = min(0.0, velocityGradient(nodeListi, i).Trace());
      const auto l = 1.0/(H(nodeListi, i).eigenValues().maxElement());
      mViscousEnergy(nodeListi, i) = (-Cl*soundSpeed(nodeListi, i) + Cq*l*mui)*l*mui;
    }
  }

  // Finally apply the boundary conditions to the viscous energy FieldList.
  for (auto boundaryItr = boundaryBegin; boundaryItr < boundaryEnd; ++boundaryItr) {
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
     const typename Dimension::Vector& /*xi*/,
     const typename Dimension::Vector& /*etai*/,
     const typename Dimension::Vector& /*vi*/,
     const typename Dimension::Scalar rhoi,
     const typename Dimension::Scalar /*csi*/,
     const typename Dimension::SymTensor& /*Hi*/,
     const typename Dimension::Vector& /*xj*/,
     const typename Dimension::Vector& /*etaj*/,
     const typename Dimension::Vector& /*vj*/,
     const typename Dimension::Scalar rhoj,
     const typename Dimension::Scalar /*csj*/,
     const typename Dimension::SymTensor& /*Hj*/) const {

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
