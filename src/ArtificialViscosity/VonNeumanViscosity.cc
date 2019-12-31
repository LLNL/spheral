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

  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  // Let the base class do it's thing.
  ArtificialViscosity<Dimension>::initialize(dataBase, state, derivs, boundaryBegin, boundaryEnd, time, dt, W);

  // Make sure the RK corrections have had boundaries completed.
  for (auto boundItr = boundaryBegin; boundItr < boundaryEnd; ++boundItr) (*boundItr)->finalizeGhostBoundary();

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
  const FieldList<Dimension, Scalar> vol = mass/massDensity;

  const RKOrder correctionOrder = this->QcorrectionOrder();

  // We'll compute the higher-accuracy RK gradient.
  FieldList<Dimension, Scalar> m0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::m0_CRKSPH);
  FieldList<Dimension, Vector> m1 = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::m1_CRKSPH);
  FieldList<Dimension, SymTensor> m2 = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::m2_CRKSPH);
  FieldList<Dimension, ThirdRankTensor> m3(FieldStorageType::CopyFields);
  FieldList<Dimension, FourthRankTensor> m4(FieldStorageType::CopyFields);
  FieldList<Dimension, Vector> gradm0 = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::gradM0_CRKSPH);
  FieldList<Dimension, Tensor> gradm1 = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::gradM1_CRKSPH);
  FieldList<Dimension, ThirdRankTensor> gradm2 = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::gradM2_CRKSPH);
  FieldList<Dimension, FourthRankTensor> gradm3(FieldStorageType::CopyFields);
  FieldList<Dimension, FifthRankTensor> gradm4(FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> A = dataBase.newFluidFieldList(0.0, "Q A");
  FieldList<Dimension, Vector> B;
  FieldList<Dimension, Tensor> C;
  FieldList<Dimension, Vector> gradA = dataBase.newFluidFieldList(Vector::zero, "Q grad A");
  FieldList<Dimension, Tensor> gradB;
  FieldList<Dimension, ThirdRankTensor> gradC;
  FieldList<Dimension, int> surfacePoint;
  if (correctionOrder == RKOrder::LinearOrder or correctionOrder == RKOrder::QuadraticOrder) {
    B = dataBase.newFluidFieldList(Vector::zero, "Q B");
    gradB = dataBase.newFluidFieldList(Tensor::zero, "Q grad B");
  }
  if (correctionOrder == RKOrder::QuadraticOrder) {
    m3 = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::m3_CRKSPH);
    m4 = dataBase.newFluidFieldList(FourthRankTensor::zero, HydroFieldNames::m4_CRKSPH);
    gradm3 = dataBase.newFluidFieldList(FourthRankTensor::zero, HydroFieldNames::gradM3_CRKSPH);
    gradm4 = dataBase.newFluidFieldList(FifthRankTensor::zero, HydroFieldNames::gradM4_CRKSPH);
    C = dataBase.newFluidFieldList(Tensor::zero, "Q C");
    gradC = dataBase.newFluidFieldList(ThirdRankTensor::zero, "Q grad C");
  }
  computeCRKSPHMoments(connectivityMap, W, vol, position, H,correctionOrder, NodeCoupling(), 
                       m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
  computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, H, surfacePoint,
                           correctionOrder, A, B, C, gradA, gradB, gradC);
  const FieldList<Dimension, Tensor> velocityGradient = gradientCRKSPH(velocity,
                                                                       position,
                                                                       vol,
                                                                       H,
                                                                       A,
                                                                       B,
                                                                       C,
                                                                       gradA,
                                                                       gradB,
                                                                       gradC,
                                                                       connectivityMap,
                                                                       correctionOrder,
                                                                       W,
                                                                       NodeCoupling());

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
