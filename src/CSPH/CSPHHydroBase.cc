//---------------------------------Spheral++----------------------------------//
// Hydro -- The CSPH/ACSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "CSPHHydroBase.hh"
#include "CSPHUtilities.hh"
#include "computeCSPHSumMassDensity.hh"
#include "computeCSPHCorrections.hh"
#include "computeHVolumes.hh"
#include "centerOfMass.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "CSPHSpecificThermalEnergyPolicy.hh"
#include "HVolumePolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/SpheralFunctions.hh"
#include "FileIO/FileIO.hh"

#include "SPH/computeSPHSumMassDensity.hh"
#include "gradientCSPH.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"

namespace Spheral {
namespace CSPHSpace {

using namespace std;
using PhysicsSpace::GenericHydro;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using Geometry::innerProduct;
using Geometry::outerProduct;

using PhysicsSpace::MassDensityType;
using PhysicsSpace::HEvolutionType;

// Put local helper methods in an anonymous namespace to keep 'em local.
namespace {

//------------------------------------------------------------------------------
// A functor for finding the equilibrium point between two kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
struct KernelFunctor {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef std::pair<double, double> ResultType;

  TableKernel<Dimension> mWT;
  double mAi, mAj, mHdeti, mHdetj, mrj, mhiinv, mhjinv, mPartialrAi, mPartialrAj;
  Vector mrjihat;

  KernelFunctor(const TableKernel<Dimension>& WT,
                const Vector& xi,
                const SymTensor& Hi,
                const Scalar Hdeti,
                const Scalar Ai,
                const Vector& gradAi,
                const Vector& xj,
                const SymTensor& Hj,
                const Scalar Hdetj,
                const Scalar Aj,
                const Vector& gradAj):
    mWT(WT),
    mAi(Ai),
    mAj(Aj),
    mHdeti(Hdeti),
    mHdetj(Hdetj),
    mrj(0.0),
    mhiinv(0.0),
    mhjinv(0.0),
    mPartialrAi(0.0),
    mPartialrAj(0.0),
    mrjihat((xj - xi).unitVector()) {
    mrj = (xj - xi).dot(mrjihat);
    mhiinv = (Hi*mrjihat).magnitude();
    mhjinv = (Hj*mrjihat).magnitude();
    mPartialrAi = gradAi.dot(mrjihat);
    mPartialrAj = gradAj.dot(mrjihat);
  }

  ResultType operator()(const double r) const {
    const double etai = mhiinv * r;
    const double etaj = mhjinv * (r - mrj);
    const std::pair<double, double> WWi = mWT.kernelAndGradValue(std::abs(etai), mHdeti);
    const std::pair<double, double> WWj = mWT.kernelAndGradValue(std::abs(etaj), mHdetj);
    const double Wi = WWi.first;
    const double Wj = WWj.first;
    const double partialrWi = mhiinv*sgn(etai)*WWi.second;
    const double partialrWj = mhjinv*sgn(etaj)*WWj.second;
    return ResultType(mAi*Wj - mAj*Wi,
                      mAi*partialrWj + mPartialrAi*Wj -
                      mAj*partialrWi - mPartialrAj*Wi);
  }
};

//------------------------------------------------------------------------------
// Search for the crossing point(s) of two kernels.
// Return value is the number of roots.  
//   Calculated root positions are x1, x2.
//   Calculated areas are dA1, dA2.
// Note: in this version we are explicitly assuming at most an A correction
// on the kernel, so B=0!
//------------------------------------------------------------------------------
template<typename Dimension>
int
kernelIntersect(const TableKernel<Dimension>& WT,
                const typename Dimension::Vector& xi,
                const typename Dimension::SymTensor& Hi,
                const typename Dimension::Scalar Hdeti,
                const typename Dimension::Scalar A0i,
                const typename Dimension::Vector& gradA0i,
                const typename Dimension::Scalar& weighti,
                const typename Dimension::Vector& xj,
                const typename Dimension::SymTensor& Hj,
                const typename Dimension::Scalar Hdetj,
                const typename Dimension::Scalar A0j,
                const typename Dimension::Vector& gradA0j,
                const typename Dimension::Scalar& weightj,
                typename Dimension::Vector& x1,
                typename Dimension::Vector& dA1,
                typename Dimension::Vector& x2,
                typename Dimension::Vector& dA2) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // // Blago!
  // x1 = 0.5*(xi + xj);
  // dA1 = (xj - xi).unitVector();
  // x2 = Vector::zero;
  // dA2 = Vector::zero;
  // return 1;
  // // Blago!

  // Build the kernel functor we use for calling into the Newton-Raphson root finder.
  KernelFunctor<Dimension> Wfunc(WT, 
                                 xi, Hi, Hdeti, A0i*weightj, gradA0i,
                                 xj, Hj, Hdetj, A0j*weighti, gradA0j);

  // We know there is a root somewhere on the j side of i.
  const double etamax = WT.kernelExtent();
  x1 = xi + newtonRaphson(Wfunc, 0.0, etamax/Wfunc.mhiinv)*Wfunc.mrjihat;
  dA1 = Wfunc.mrjihat;   // Currently hard-wired for 1D!

  // Check for a root on the backside of i.
  if (etamax/Wfunc.mhiinv < etamax/Wfunc.mhjinv - Wfunc.mrj) {

    x2 = xi + newtonRaphson(Wfunc, -etamax/Wfunc.mhiinv, 0.0)*Wfunc.mrjihat;
    dA2 = -Wfunc.mrjihat;
    return 2;

  } else if (etamax/Wfunc.mhjinv < etamax/Wfunc.mhiinv - Wfunc.mrj) {

    // Similarly check for a root on the backside of j.
    x2 = xi + newtonRaphson(Wfunc, Wfunc.mrj, Wfunc.mrj + etamax/Wfunc.mhjinv)*Wfunc.mrjihat;
    dA2 = -Wfunc.mrjihat;
    return 2;

  }

  // There was only one root.
  x2 = Vector::zero;
  dA2 = Vector::zero;
  return 1;
}

//------------------------------------------------------------------------------
// Compute the net pair-wise force.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
pairWiseForce(const TableKernel<Dimension>& WT,
              const typename Dimension::Vector& xi,
              const typename Dimension::SymTensor& Hi,
              const typename Dimension::Scalar Hdeti,
              const typename Dimension::Scalar A0i,
              const typename Dimension::Vector& gradA0i,
              const typename Dimension::Scalar& weighti,
              const typename Dimension::Scalar& Pi,
              const typename Dimension::Tensor& Qi,
              const typename Dimension::Vector& xj,
              const typename Dimension::SymTensor& Hj,
              const typename Dimension::Scalar Hdetj,
              const typename Dimension::Scalar A0j,
              const typename Dimension::Vector& gradA0j,
              const typename Dimension::Scalar& weightj,
              const typename Dimension::Scalar& Pj,
              const typename Dimension::Tensor& Qj) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Find the point(s) and area(s) where the kernels of the two points are equal.
  // This defines the surface for the volumes of the two points.
  Vector x1, x2, dA1, dA2;
  const int nsurf = kernelIntersect(WT, 
                                    xi, Hi, Hdeti, A0i, gradA0i, weighti,
                                    xj, Hj, Hdetj, A0j, gradA0j, weightj,
                                    x1, dA1, x2, dA2);
  CHECK(nsurf == 1 or nsurf == 2);

  // Weight at each point.
  const Scalar wj = weightj*A0i*WT.kernelValue((Hj*(xi - xj)).magnitude(), Hdetj);
  const Scalar wi = weighti*A0j*WT.kernelValue((Hi*(xj - xi)).magnitude(), Hdeti);

  // // Determine the etas, gradP, and gradQ.
  // const Vector etai = Hi*(xj - xi);
  // const Vector Pgrad = (Pj - Pi)*etai/std::max(1.0e-30, etai.magnitude2());
  // const ThirdRankTensor Qgrad = outerProduct<Dimension>(Qj - Qi, etai/std::max(1.0e-30, etai.magnitude2()));

  // // Linearly interpolate the pressure to the intersection points, and sum
  // // the force.
  // const Scalar P1 = Pi + Pgrad.dot(Hi*(x1 - xi));
  // const Tensor Q1 = Qi + innerProduct<Dimension>(Qgrad, Hi*(x1 - xi));
  // // const Scalar wj1 = weightj*A0i*WT.kernelValue((Hj*(x1 - xj)).magnitude(), Hdetj);
  // // const Scalar wi1 = weighti*A0j*WT.kernelValue((Hi*(x1 - xi)).magnitude(), Hdeti);
  // const Scalar wij1 = 0.5*(wj + wi);
  // Vector result = -wij1*(P1*dA1 + Q1*dA1);

  // // Is there a second intersection?
  // if (nsurf == 2) {
  //   const Scalar P2 = Pi + Pgrad.dot(Hi*(x2 - xi));
  //   const Tensor Q2 = Qi + innerProduct<Dimension>(Qgrad, Hi*(x2 - xi));
  //   // const Scalar wj2 = weightj*A0i*WT.kernelValue((Hj*(x2 - xj)).magnitude(), Hdetj);
  //   // const Scalar wi2 = weighti*A0j*WT.kernelValue((Hi*(x2 - xi)).magnitude(), Hdeti);
  //   const Scalar wij2 = 0.5*(wj + wi);
  //   result -= wij2*(P2*dA2 + Q2*dA2);
  // }

  // Sum the pressure at the effective face.
  const Scalar wj1 = A0i*weightj*WT.kernelValue((Hj*(xi - xj)).magnitude(), Hdetj);
  const Scalar wi1 = A0j*weighti*WT.kernelValue((Hi*(xj - xi)).magnitude(), Hdeti);
  Vector result = -(wi1*Pi + wj1*Pj)*dA1 - (wi1*Qi + wj1*Qj)*dA1;

  // Is there a second intersection?
  if (nsurf == 2) {
    const Scalar wj2 = A0i*weightj*WT.kernelValue((Hj*(xi - xj)).magnitude(), Hdetj);
    const Scalar wi2 = A0j*weighti*WT.kernelValue((Hi*(xj - xi)).magnitude(), Hdeti);
    result -= (wi2*Pi + wj2*Pj)*dA2 + (wi2*Qi + wj2*Qj)*dA2;
  }

  // That's it.
  return result;
}

}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
CSPHHydroBase<Dimension>::
CSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
              const TableKernel<Dimension>& W,
              const TableKernel<Dimension>& WPi,
              ArtificialViscosity<Dimension>& Q,
              const double filter,
              const double cfl,
              const bool useVelocityMagnitudeForDt,
              const bool compatibleEnergyEvolution,
              const bool XSPH,
              const MassDensityType densityUpdate,
              const HEvolutionType HUpdate):
  GenericHydro<Dimension>(W, WPi, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mXSPH(XSPH),
  mfilter(filter),
  mTimeStepMask(FieldSpace::Copy),
  mPressure(FieldSpace::Copy),
  mSoundSpeed(FieldSpace::Copy),
  mVolume(FieldSpace::Copy),
  mSpecificThermalEnergy0(FieldSpace::Copy),
  mHideal(FieldSpace::Copy),
  mMaxViscousPressure(FieldSpace::Copy),
  // mMassDensitySum(FieldSpace::Copy),
  mWeightedNeighborSum(FieldSpace::Copy),
  mMassSecondMoment(FieldSpace::Copy),
  mXSPHDeltaV(FieldSpace::Copy),
  mDxDt(FieldSpace::Copy),
  mDvDt(FieldSpace::Copy),
  mDmassDensityDt(FieldSpace::Copy),
  mDspecificThermalEnergyDt(FieldSpace::Copy),
  mDHDt(FieldSpace::Copy),
  mDvDx(FieldSpace::Copy),
  mInternalDvDx(FieldSpace::Copy),
  mDmassDensityDx(FieldSpace::Copy),
  mPairAccelerations(FieldSpace::Copy),
  mM0(FieldSpace::Copy),
  mM1(FieldSpace::Copy),
  mM2(FieldSpace::Copy),
  mA0(FieldSpace::Copy),
  mA(FieldSpace::Copy),
  mB(FieldSpace::Copy),
  mC(FieldSpace::Copy),
  mD(FieldSpace::Copy),
  mGradA0(FieldSpace::Copy),
  mGradA(FieldSpace::Copy),
  mGradB(FieldSpace::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
CSPHHydroBase<Dimension>::
~CSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
  // mMassDensitySum = dataBase.newFluidFieldList(0.0, ReplaceFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity);
  mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mDmassDensityDx = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::massDensityGradient);
  mPairAccelerations = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::pairAccelerations);

  mM0 = dataBase.newFluidFieldList(0.0,             HydroFieldNames::m0_CSPH);
  mM1 = dataBase.newFluidFieldList(Vector::zero,    HydroFieldNames::m1_CSPH);
  mM2 = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::m2_CSPH);
  mA0 = dataBase.newFluidFieldList(0.0,             HydroFieldNames::A0_CSPH);
  mA = dataBase.newFluidFieldList(0.0,              HydroFieldNames::A_CSPH);
  mB = dataBase.newFluidFieldList(Vector::zero,     HydroFieldNames::B_CSPH);
  mC = dataBase.newFluidFieldList(Vector::zero,     HydroFieldNames::C_CSPH);
  mD = dataBase.newFluidFieldList(Tensor::zero,     HydroFieldNames::D_CSPH);
  mGradA0 = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::gradA0_CSPH);
  mGradA = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::gradA_CSPH);
  mGradB = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::gradB_CSPH);

  // // Compute the volumes.
  // const TableKernel<Dimension>& W = this->kernel();
  // const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  // computeHVolumes(W.kernelExtent(), H, mVolume);

  // // We need the boundary information on volume to initialize the CSPH corrections.
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->applyFieldListGhostBoundary(mVolume);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // // Compute the kernel correction fields.
  // const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  // computeCSPHCorrections(connectivityMap, W, mVolume, position, H, mM0, mM1, mM2, mA0, mA, mB, mC, mD, mGradA, mGradB);

  // Initialize the pressure and sound speed.
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Create the local storage for time step mask, pressure, sound speed, and correction fields.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  // dataBase.fluidPressure(mPressure);
  // dataBase.fluidSoundSpeed(mSoundSpeed);
  dataBase.resizeFluidFieldList(mPressure, 0.0,          HydroFieldNames::pressure, false);
  dataBase.resizeFluidFieldList(mSoundSpeed, 0.0,        HydroFieldNames::soundSpeed, false);
  dataBase.resizeFluidFieldList(mM0,    0.0,             HydroFieldNames::m0_CSPH, false);
  dataBase.resizeFluidFieldList(mM1,    Vector::zero,    HydroFieldNames::m1_CSPH, false);
  dataBase.resizeFluidFieldList(mM2,    SymTensor::zero, HydroFieldNames::m2_CSPH, false);
  dataBase.resizeFluidFieldList(mA0,    0.0,             HydroFieldNames::A0_CSPH, false);
  dataBase.resizeFluidFieldList(mA,     0.0,             HydroFieldNames::A_CSPH, false);
  dataBase.resizeFluidFieldList(mB,     Vector::zero,    HydroFieldNames::B_CSPH, false);
  dataBase.resizeFluidFieldList(mC,     Vector::zero,    HydroFieldNames::C_CSPH, false);
  dataBase.resizeFluidFieldList(mD,     Tensor::zero,    HydroFieldNames::D_CSPH, false);
  dataBase.resizeFluidFieldList(mGradA0,Vector::zero,    HydroFieldNames::gradA0_CSPH, false);
  dataBase.resizeFluidFieldList(mGradA, Vector::zero,    HydroFieldNames::gradA_CSPH, false);
  dataBase.resizeFluidFieldList(mGradB, Tensor::zero,    HydroFieldNames::gradB_CSPH, false);

  // If we're using the compatibile energy discretization, prepare to maintain a copy
  // of the thermal energy.
  dataBase.resizeFluidFieldList(mSpecificThermalEnergy0, 0.0);
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      *mSpecificThermalEnergy0[nodeListi] = (*itr)->specificThermalEnergy();
      (*mSpecificThermalEnergy0[nodeListi]).name(HydroFieldNames::specificThermalEnergy + "0");
    }
  }

  // Now register away.
  // Mass.
  FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  state.enroll(mass);

  // Volume.
  const TableKernel<Dimension>& W = this->kernel();
  // PolicyPointer volumePolicy(new HVolumePolicy<Dimension>(W.kernelExtent()));
  // state.enroll(mVolume, volumePolicy);

  // We need to build up CompositeFieldListPolicies for the mass density and H fields
  // in order to enforce NodeList dependent limits.
  FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  boost::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > rhoPolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  boost::shared_ptr<CompositeFieldListPolicy<Dimension, SymTensor> > Hpolicy(new CompositeFieldListPolicy<Dimension, SymTensor>());
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    rhoPolicy->push_back(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                      (*itr)->rhoMax()));
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == PhysicsSpace::IntegrateH) {
      Hpolicy->push_back(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    } else {
      CHECK(HEvolution() == PhysicsSpace::IdealH);
      Hpolicy->push_back(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    }
  }
  state.enroll(massDensity, rhoPolicy);
  state.enroll(Hfield, Hpolicy);

  // Register the position update, which depends on whether we're using XSPH or not.
  FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  if (mXSPH) {
    PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
    state.enroll(position, positionPolicy);
  } else {
    PolicyPointer positionPolicy(new PositionPolicy<Dimension>());
    state.enroll(position, positionPolicy);
  }

  // Are we using the compatible energy evolution scheme?
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  if (compatibleEnergyEvolution()) {
    PolicyPointer thermalEnergyPolicy(new SpecificThermalEnergyPolicy<Dimension>(dataBase));
    // PolicyPointer thermalEnergyPolicy(new NonSymmetricSpecificThermalEnergyPolicy<Dimension>(dataBase));
    // PolicyPointer thermalEnergyPolicy(new CSPHSpecificThermalEnergyPolicy<Dimension>(dataBase, this->kernel()));
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,
                                                                           HydroFieldNames::specificThermalEnergy));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
    state.enroll(mSpecificThermalEnergy0);
  } else {
    PolicyPointer thermalEnergyPolicy(new IncrementFieldList<Dimension, Scalar>());
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>());
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
  }

  // Register the time step mask, initialized to 1 so that everything defaults to being
  // checked.
  state.enroll(mTimeStepMask);

  // Compute and register the pressure and sound speed.
  PolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
  PolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mSoundSpeed, csPolicy);

  // Register the CSPH correction fields.
  // We deliberately make these non-dynamic here.  This corrections are computed
  // during CSPHHydroBase::initialize, not as part of our usual state update.
  // This is necessary 'cause we need boundary conditions *and* the current set of
  // neighbors before we compute these suckers.
  state.enroll(mM0);
  state.enroll(mM1);
  state.enroll(mM2);
  state.enroll(mA0);
  state.enroll(mA);
  state.enroll(mB);
  state.enroll(mC);
  state.enroll(mD);
  state.enroll(mGradA0);
  state.enroll(mGradA);
  state.enroll(mGradB);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  typedef typename StateDerivatives<Dimension>::KeyType Key;
  const string DxDtName = IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const string DvDtName = IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  // dataBase.resizeFluidFieldList(mMassDensitySum, 0.0, ReplaceFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mDmassDensityDx, Vector::zero, HydroFieldNames::massDensityGradient, false);
  dataBase.resizeFluidFieldList(mPairAccelerations, vector<Vector>(), HydroFieldNames::pairAccelerations, false);

  derivs.enroll(mHideal);
  derivs.enroll(mMaxViscousPressure);
  // derivs.enroll(mMassDensitySum);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mXSPHDeltaV);

  // These two (the position and velocity updates) may be registered
  // by other physics packages as well, so we need to be careful
  // not to duplicate if so.
  if (not derivs.registered(mDxDt)) derivs.enroll(mDxDt);
  if (not derivs.registered(mDvDt)) derivs.enroll(mDvDt);

  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDHDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enroll(mDmassDensityDx);
  derivs.enroll(mPairAccelerations);
}

//------------------------------------------------------------------------------
// Initialize the hydro before evaluating derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // Get the artificial viscosity and initialize it.
  const TableKernel<Dimension>& W = this->kernel();
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();
  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               W);

  // Compute the kernel correction fields.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CSPH, 0.0);
  FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CSPH, Vector::zero);
  FieldList<Dimension, SymTensor> m2 = state.fields(HydroFieldNames::m2_CSPH, SymTensor::zero);
  FieldList<Dimension, Scalar> A0 = state.fields(HydroFieldNames::A0_CSPH, 0.0);
  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  FieldList<Dimension, Vector> gradA0 = state.fields(HydroFieldNames::gradA0_CSPH, Vector::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);

  // Change CSPH weights here if need be!
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> vol = mass/massDensity;
  computeCSPHCorrections(connectivityMap, W, vol, position, H, true, m0, m1, m2, A0, A, B, C, D, gradA0, gradA, gradB);
  // computeCSPHCorrections(connectivityMap, W, mass, position, H, true, m0, m1, m2, A0, A, B, C, D, gradA0, gradA, gradB);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(A0);
    (*boundItr)->applyFieldListGhostBoundary(A);
    (*boundItr)->applyFieldListGhostBoundary(B);
    (*boundItr)->applyFieldListGhostBoundary(C);
    (*boundItr)->applyFieldListGhostBoundary(D);
    (*boundItr)->applyFieldListGhostBoundary(gradA0);
    (*boundItr)->applyFieldListGhostBoundary(gradA);
    (*boundItr)->applyFieldListGhostBoundary(gradB);
  }

  // // If we're doing the RigorousSumDensity update, now is a good time to do it
  // // since we have the boundary conditions and corrections all ready to go.
  // if (densityUpdate() == PhysicsSpace::RigorousSumDensity) {
  //   const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  //   FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  //   computeCSPHSumMassDensity(connectivityMap, W, position, mass, vol, H, A0, A, B, massDensity);
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
  // }

}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  // const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CSPH, 0.0);
  const FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CSPH, Vector::zero);
  const FieldList<Dimension, Scalar> A0 = state.fields(HydroFieldNames::A0_CSPH, 0.0);
  const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  const FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  const FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  const FieldList<Dimension, Vector> gradA0 = state.fields(HydroFieldNames::gradA0_CSPH, Vector::zero);
  const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);

  // CHECK(vol.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(m0.size() == numNodeLists);
  CHECK(m1.size() == numNodeLists);
  CHECK(A0.size() == numNodeLists);
  CHECK(A.size() == numNodeLists);
  CHECK(B.size() == numNodeLists);
  CHECK(C.size() == numNodeLists);
  CHECK(D.size() == numNodeLists);
  CHECK(gradA0.size() == numNodeLists);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists);

  // // BLAGO!  As a sanity check recompute all the CSPH corrections now.
  // computeCSPHCorrections(connectivityMap, W, vol, position, H, m0, m1, const_cast<FieldList<Dimension, SymTensor>&>(mM2), A0, A, B, C, D, gradA, gradB);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) {
  //   (*boundItr)->applyFieldListGhostBoundary(A0);
  //   (*boundItr)->applyFieldListGhostBoundary(A);
  //   (*boundItr)->applyFieldListGhostBoundary(B);
  //   (*boundItr)->applyFieldListGhostBoundary(C);
  //   (*boundItr)->applyFieldListGhostBoundary(D);
  //   (*boundItr)->applyFieldListGhostBoundary(gradA);
  //   (*boundItr)->applyFieldListGhostBoundary(gradB);
  // }
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // // BLAGO!

  // Derivative FieldLists.
  // FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, Vector> DrhoDx = derivatives.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  // CHECK(rhoSum.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(DrhoDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const NodeList<Dimension>& nodeList = **itr;
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar hmin = nodeList.hmin();
    const Scalar hmax = nodeList.hmax();
    const Scalar hminratio = nodeList.hminratio();
    const int maxNumNeighbors = nodeList.maxNumNeighbors();
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();

    // Get the work field for this NodeList.
    Field<Dimension, Scalar>& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Scalar& epsi = specificThermalEnergy(nodeListi, i);
      const Scalar& Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar& ci = soundSpeed(nodeListi, i);
      const Scalar& m0i = m0(nodeListi, i);
      const Vector& m1i = m1(nodeListi, i);
      const Scalar& A0i = A0(nodeListi, i);
      const Scalar& Ai = A(nodeListi, i);
      const Vector& Bi = B(nodeListi, i);
      const Vector& gradA0i = gradA0(nodeListi, i);
      const Vector& gradAi = gradA(nodeListi, i);
      const Tensor& gradBi = gradB(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar weighti = mi/rhoi;  // Change CSPH weights here if need be!
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Ai > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(weighti > 0.0);

      // Scalar& rhoSumi = rhoSum(nodeListi, i);
      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      Vector& DrhoDxi = DrhoDx(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Bizarrely, in CSPH there is a self-contribution to gradients.  We need this 
      // term to compute those.
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      const Vector selfGradContrib = W0*(Ai*Bi + gradAi);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = 1.0; // (nodeListi == nodeListj ? 1.0 : 0.2);
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const Vector& rj = position(nodeListj, j);
              const Scalar& mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar& rhoj = massDensity(nodeListj, j);
              const Scalar& epsj = specificThermalEnergy(nodeListj, j);
              const Scalar& Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar& cj = soundSpeed(nodeListj, j);
              const Scalar& A0j = A0(nodeListj, j);
              const Scalar& Aj = A(nodeListj, j);
              const Vector& Bj = B(nodeListj, j);
              const Vector& gradA0j = gradA0(nodeListj, j);
              const Vector& gradAj = gradA(nodeListj, j);
              const Tensor& gradBj = gradB(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar weightj = mj/rhoj;     // Change CSPH weights here if need be!
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Aj > 0.0 or j >= firstGhostNodej);
              CHECK(Hdetj > 0.0);
              CHECK(weightj > 0.0);

              // Scalar& rhoSumj = rhoSum(nodeListj, j);
              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Vector& DrhoDxj = DrhoDx(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              Vector& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              Scalar& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              SymTensor& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Node displacement.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);

              // if (i == 0) {
              //   cerr << j << " " << rj << " " << Pj << " " << -rij << " " << -etaj << endl;
              // }

              // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);

              // Symmetrized kernel weight and gradient.
              Scalar gWi, gWj, W1i, W1j, W0i, W0j;
              Vector gradW1i, gradW1j, gradW0i, gradW0j;
              CSPHKernelAndGradient(W,  rij,  etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, W1j, gWj, gradW1j);
              CSPHKernelAndGradient(W, -rij, -etai, Hi, Hdeti, Aj, Bj, gradAj, gradBj, W1i, gWi, gradW1i);
              CSPHKernelAndGradient(W,  rij,  etaj, Hj, Hdetj, A0i, Vector::zero, gradA0i, Tensor::zero, W0j, gWj, gradW0j);
              CSPHKernelAndGradient(W, -rij, -etai, Hi, Hdeti, A0j, Vector::zero, gradA0j, Tensor::zero, W0i, gWi, gradW0i);
              const Vector gradWSPHi = gWi*(Hi*etai.unitVector());
              const Vector gradWSPHj = gWj*(Hj*etaj.unitVector());

              // Compute the limiter determining how much of the linear correction we allow.
              // const Scalar fQ = max(0.0, min(1.0, min(max(0.0, abs(0.05*Pi) - Qi)*safeInv(abs(0.05*Pi)), 
              //                                         max(0.0, abs(0.05*Pj) - Qj)*safeInv(abs(0.05*Pj)))));
              // const Scalar fL = max(0.0, min(1.0, 1.0 - abs(0.5*(DrhoDxj + DrhoDxi).dot(rij) - (rhoi - rhoj))));
              const Scalar fg = max(0.0, min(1.0, -(gradW1j.dot(gradW1i))*safeInv(sqrt(gradW1j.magnitude2()*gradW1i.magnitude2()))));
              const Scalar f = fg; //min(fQ, min(fL, fg));
              CHECK2(f >= 0.0 and f <= 1.0, "Failing f bounds: " << f);
              const Scalar Wj = (1.0 - f)*W0j + f*W1j;
              const Scalar Wi = (1.0 - f)*W0i + f*W1i;
              const Vector gradWj = (1.0 - f)*gradW0j + f*gradW1j;
              const Vector gradWi = (1.0 - f)*gradW0i + f*gradW1i;
              // CHECK(gradWj.dot(gradWi) <= 0.0);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()/(rij2 + 1.0e-10) / FastMath::square(Dimension::pownu12(rij2 + 1.0e-10));
              weightedNeighborSumi += fweightij*std::abs(gWi);
              weightedNeighborSumj += fweightij*std::abs(gWj);
              massSecondMomenti += fweightij*gradWSPHi.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWSPHj.magnitude2()*thpt;

              // Velocity gradient.
              const Vector vij = vi - vj;
              const Tensor deltaDvDxi = weightj*vj.dyad(gradWj);
              const Tensor deltaDvDxj = weighti*vi.dyad(gradWi);
              DvDxi += deltaDvDxi;
              DvDxj += deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi += deltaDvDxi;
                localDvDxj += deltaDvDxj;
              }

              // Mass density gradient.
              DrhoDxi += weightj*rhoj*gradWj;
              DrhoDxj += weighti*rhoi*gradWi;

              // Acceleration (CSPH form).
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const Vector deltagrad = gradWj - gradWi;
              // const Vector forceij = 0.5*weighti*weightj*((Pi + Pj)*deltagrad + 
              //                                             (rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second)*deltagrad);
              const Vector forceij = 0.5*weighti*weightj*(Pi + Pj)*deltagrad + 0.5*(mi*mj*QPiij.first*gradWSPHi + mi*mj*QPiij.second*gradWSPHj);
              Vector deltaDvDti = -forceij/mi;
              Vector deltaDvDtj =  forceij/mj;
              DvDti += deltaDvDti;
              DvDtj += deltaDvDtj;
              if (mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(deltaDvDti);
                pairAccelerationsj.push_back(deltaDvDtj);
              }

              // Specific thermal energy evolution.
              // Equipartition the q work.
              // const Scalar workQ = (rhoi*rhoi*(QPiij.first*vij).dot(deltagrad0) +
              //                       rhoj*rhoj*(QPiij.second*vij).dot(deltagrad0));
              // DepsDti += 0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + 0.5*workQ)/mi;
              // DepsDtj += 0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + 0.5*workQ)/mj;

              // Q work based on the Q per point.
              // DepsDti += 0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + rhoi*rhoi*(QPiij.first*vij).dot(deltagrad0))/mi;
              // DepsDtj += 0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + rhoj*rhoj*(QPiij.second*vij).dot(deltagrad0))/mj;
              // DepsDti += 0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + rhoj*rhoj*(QPiij.second*vij).dot(deltagrad0))/mi;
              // DepsDtj += 0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + rhoi*rhoi*(QPiij.first*vij).dot(deltagrad0))/mj;

              // CSPH PdV and SPH Q work.
              DepsDti += 0.5*weighti*weightj*Pj*vij.dot(deltagrad)/mi + 0.25*mj*(QPiij.first *vij).dot(gradWSPHi + gradWSPHj);
              DepsDtj += 0.5*weighti*weightj*Pi*vij.dot(deltagrad)/mj + 0.25*mi*(QPiij.second*vij).dot(gradWSPHi + gradWSPHj);
              
              // Estimate of delta v (for XSPH).
              if (mXSPH and (nodeListi == nodeListj)) {
                XSPHDeltaVi -= weightj*Wj*vij;
		XSPHDeltaVj += weighti*Wi*vij;
              }

            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not mCompatibleEnergyEvolution or 
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // // Finish the mass density sum.
      // rhoSumi = A0i*(rhoSumi + mi*W(0.0, Hdeti));
      // // rhoSumi += mi*W(0.0, Hdeti);

      // Finish the velocity gradient.
      DvDxi += weighti*vi*selfGradContrib;
      localDvDxi += weighti*vi*selfGradContrib;

      // Finish the density gradient.
      DrhoDxi += weighti*rhoi*selfGradContrib;

      // Time evolution of the mass density.
      DrhoDti = -rhoi*DvDxi.Trace();

      // Finish the acceleration.
      // const Vector deltaDvDtii = -weighti/rhoi*Pi*selfGradContrib;
      // const Vector deltaDvDtii = weighti*Pi*Ai*W0*Bi.unitVector();
      // DvDti += deltaDvDtii;

      // // The specific thermal energy evolution.
      // DepsDti -= Pi/rhoi*DvDxi.Trace();

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             ri,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                        ri,
                                                        weightedNeighborSumi,
                                                        massSecondMomenti,
                                                        W,
                                                        hmin,
                                                        hmax,
                                                        hminratio,
                                                        nPerh,
                                                        connectivityMap,
                                                        nodeListi,
                                                        i);

      //const Scalar mag0 = DxDti.magnitude();
      const Scalar mag0 = vi.magnitude();
      //printf("MAG0=%10.3e",mag0);
      
      // if (dt > 0.0) {
      //   CHECK(m0i > 0.0);
      //   const Vector com = -m1i/m0i;
      //   const Vector dhat = com.unitVector();
      //   //const Vector delPos=com - ri;
      //   const Scalar a0 = DvDti.magnitude();
      //   const Vector delPos=com;
      //   //const Vector accel=2*delPos/(dt*dt)-2*vi/dt;
      //   //const Vector accel=2*delPos/(dt*dt)-2*DxDti/dt;
      //   const Scalar deltamag = com.magnitude();
      //   //const Vector accel=std::min(0.01*mag0, deltamag)*2*delPos/(dt*dt);
      //   const Vector accel=2*delPos/(dt*dt);
      //   const Scalar a1 = accel.magnitude();
      //   const Vector delta = mfilter*std::min(a0, a1)*accel.unitVector();
      //   //const Vector delta = 0.01*std::min(a0, a1)*accel.unitVector();

      //   // const Vector delPos2=std::min(0.01*mag0, deltamag)*dhat;
      //   //const Vector accel=2*delPos2/(dt*dt)-2*vi/dt;
      //   //const Vector accel=2*delPos/(dt*dt*mi);
      //   // printf("DVDT=%10.3e, accell=%10.3e, delta=%10.3e\n",DvDti[0],accel[0],delta[0]);
      //   // printf("COM=%10.3e, ri=%10.3e, del=%10.3e dt=%10.3e vi=%10.3e\n",com[0],ri[0],delPos[0],dt,vi[0]);
      //   //DvDti += accel;
      //   DvDti += delta;

      //   // Account for the work done as well.
      //   // DepsDti -= (vi + 0.5*dt*DvDti).dot(DvDti);
      // }

      // if (dt > 0.0) {
      //    const Vector com = centerOfMass(polyvol(nodeListi, i), DrhoDxi);
      //    const Vector dhat = com.unitVector();
      //    //const Vector delPos=com - ri;
      //    const Scalar a0 = DvDti.magnitude();
      //    const Vector delPos=com;
      //    //const Vector accel=2*delPos/(dt*dt)-2*vi/dt;
      //    //const Vector accel=2*delPos/(dt*dt)-2*DxDti/dt;
      //    const Scalar deltamag = com.magnitude();
      //    //const Vector accel=std::min(0.01*mag0, deltamag)*2*delPos/(dt*dt);
      //    const Vector accel=2*delPos/(dt*dt);
      //    const Scalar a1 = accel.magnitude();
      //    const Vector delta = mfilter*std::min(a0, a1)*accel.unitVector();
      //    //const Vector delta = 0.01*std::min(a0, a1)*accel.unitVector();

      //    // const Vector delPos2=std::min(0.01*mag0, deltamag)*dhat;
      //    //const Vector accel=2*delPos2/(dt*dt)-2*vi/dt;
      //    //const Vector accel=2*delPos/(dt*dt*mi);
      //    // printf("DVDT=%10.3e, accell=%10.3e, delta=%10.3e\n",DvDti[0],accel[0],delta[0]);
      //    // printf("COM=%10.3e, ri=%10.3e, del=%10.3e dt=%10.3e vi=%10.3e\n",com[0],ri[0],delPos[0],dt,vi[0]);
      //    //DvDti += accel;
      //    DvDti += delta;
      
      // }

      // // BLAGO!
      // if (i == 0 or i == 50) {
      //   cerr.precision(13);
      //   cerr << " --> " << i;
      //   if (i == 0) {
      //     cerr << "  ";
      //   } else {
      //     cerr << " ";
      //   }
      //   cerr << " " << Ai
      //        << " " << Bi
      //        << " " << rhoi
      //        // << " " << rhoSumi 
      //        << " " << DxDti 
      //        << " " << DrhoDti 
      //        << " " << DvDti
      //        << " " << DepsDti 
      //        << " " << DvDxi 
      //        << " " << DrhoDxi
      //        << " " << DHDti
      //        << " " << Hideali
      //        << " " << maxViscousPressurei
      //        << endl;
      // }
      // // BLAGO!

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());

      // Now add the pairwise time for each neighbor we computed here.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
          Field<Dimension, Scalar>& workFieldj = nodeLists[nodeListj]->work();
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              workFieldj(j) += deltaTimePair;
            }
          }
        }
      }
    }
  }

  // BLAGO!
  // {
  //   FieldList<Dimension, Tensor> DvDx_check = gradientCSPH(velocity, position, vol, H, A, B, C, D, gradA, gradB, connectivityMap, W);
  //   FieldList<Dimension, Vector> gradP_check = gradientCSPH(pressure, position, vol, H, A, B, C, D, gradA, gradB, connectivityMap, W);
  //   for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     for (size_t i = 0; i != DvDx[nodeListi]->numInternalElements(); ++i) {
  // 	if (std::abs(DvDx(nodeListi, i).xx() - DvDx_check(nodeListi, i).xx()) > 1e-6)
  // 	  cerr << "DVDX fail: (" << i << " " << DvDx(nodeListi, i) << " " << DvDx_check(nodeListi, i) << ") " << endl;
  // 	const Vector DvDti_check = -gradP_check(nodeListi, i)/massDensity(nodeListi, i);
  // 	if (std::abs(DvDt(nodeListi, i).x() - DvDti_check.x()) > 1.0e-6)
  // 	  cerr << "DVDT fail: (" << i << " " << DvDt(nodeListi, i) << " " << DvDti_check << ") " << endl;
  //     }
  //   }
  // }
  // BLAGO!

}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    FieldList<Dimension, Vector> accelerations = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Finalize the state after state has been updated and boundary conditions 
// enforced.  For CSPH this is where we update the volumes and RPKM corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
postStateUpdate(const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                const StateDerivatives<Dimension>& derivs) const {

  // // Grab state we're going to use.
  // const TableKernel<Dimension>& W = this->kernel();
  // const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  // const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  // const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);

  // // Compute the volume per node.
  // FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  // computeHullVolumes(connectivityMap, position, vol);

  // // We need boundary conditions enforced on the volume before we can compute corrections.
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->applyFieldListGhostBoundary(vol);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // // Compute the kernel correction fields.
  // FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  // FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  // FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  // FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  // FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  // FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);
  // computeCSPHCorrections(connectivityMap, W, vol, position, H, A, B, C, D, gradA, gradB);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) {
  //   (*boundItr)->applyFieldListGhostBoundary(A);
  //   (*boundItr)->applyFieldListGhostBoundary(B);
  //   (*boundItr)->applyFieldListGhostBoundary(C);
  //   (*boundItr)->applyFieldListGhostBoundary(D);
  //   (*boundItr)->applyFieldListGhostBoundary(gradA);
  //   (*boundItr)->applyFieldListGhostBoundary(gradB);
  // }
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Base class finalization.
  GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Add any filtering component to the node movement.
  // Note that the FacetedVolumes are in coordinates with the node at the origin already!
  if (mfilter > 0.0) {
    const TableKernel<Dimension>& W = this->kernel();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const FieldList<Dimension, Vector> DrhoDx = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero);
    const unsigned numNodeLists = mass.size();
    const Scalar W0 = W.kernelValue(0.0, 1.0);
    FieldList<Dimension, Vector> delta = dataBase.newFluidFieldList(Vector::zero, "delta position");
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
           iItr != connectivityMap.end(nodeListi);
           ++iItr) {
        const int i = *iItr;
        const Vector& ri = position(nodeListi, i);
        const Scalar mi = mass(nodeListi, i);
        const Scalar rhoi = massDensity(nodeListi, i);
        const Vector DrhoDxi = DrhoDx(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          for (typename vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
               jItr != fullConnectivity[nodeListj].end();
               ++jItr) {
            const unsigned j = *jItr;
            const Vector& rj = position(nodeListj, j);
            const Scalar mj = mass(nodeListj, j);
            const Scalar rhoj = massDensity(nodeListj, j);
            const Vector DrhoDxj = DrhoDx(nodeListj, j);
            const Vector rji = rj - ri;
            const Vector rjihat = rji.unitVector();
            const Scalar rhoij = rhoi + 0.25*DrhoDxi.dot(rji);
            const Scalar rhoji = rhoj - 0.25*DrhoDxj.dot(rji);
            const Scalar deltaj = max(0.0, 0.5*(Dimension::rootnu(mi/rhoij) + Dimension::rootnu(mj/rhoji)) - rji.magnitude());
            const Scalar etai = (Hi*rji).magnitude();
            const Scalar weight = W.kernelValue(etai, 1.0)/W0;
            delta(nodeListi, i) -= weight*deltaj*rjihat;
          }
        }
      }
    }

    // Apply the filtering.
    const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
        if (mag0 > 0.0) {
          const Scalar deltamag = delta(nodeListi, i).magnitude();
          const Scalar effmag = min(mfilter*mag0, deltamag);
          position(nodeListi, i) += effmag*delta(nodeListi, i).unitVector();
        }
      }
    }

    // Check for any boundary violations.
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
    this->enforceBoundaries(state, derivs);
  }

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (densityUpdate() == PhysicsSpace::RigorousSumDensity) {
    const TableKernel<Dimension>& W = this->kernel();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    // FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);

    const FieldList<Dimension, Scalar> vol = mass/massDensity;

    computeCSPHSumMassDensity(connectivityMap, this->kernel(), position, mass, vol, H, this->boundaryBegin(), this->boundaryEnd(), massDensity);
    // SPHSpace::computeSPHSumMassDensity(connectivityMap, this->kernel(), position, mass, H, massDensity);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // } else if (densityUpdate() == PhysicsSpace::SumDensity) {
  //   FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   FieldList<Dimension, Scalar> massDensitySum = derivs.fields(ReplaceFieldList<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
  //                                                               HydroFieldNames::massDensity, 0.0);
  //   massDensity.assignFields(massDensitySum);
  }

}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Apply boundary conditions to the basic fluid state Fields.
  // FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);

  FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CSPH, 0.0);
  FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CSPH, Vector::zero);
  FieldList<Dimension, Scalar> A0 = state.fields(HydroFieldNames::A0_CSPH, 0.0);
  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  FieldList<Dimension, Vector> gradA0 = state.fields(HydroFieldNames::gradA0_CSPH, Vector::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    // (*boundaryItr)->applyFieldListGhostBoundary(vol);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    if (compatibleEnergyEvolution()) (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
    (*boundaryItr)->applyFieldListGhostBoundary(m0);
    (*boundaryItr)->applyFieldListGhostBoundary(m1);
    (*boundaryItr)->applyFieldListGhostBoundary(A0);
    (*boundaryItr)->applyFieldListGhostBoundary(A);
    (*boundaryItr)->applyFieldListGhostBoundary(B);
    (*boundaryItr)->applyFieldListGhostBoundary(C);
    (*boundaryItr)->applyFieldListGhostBoundary(D);
    (*boundaryItr)->applyFieldListGhostBoundary(gradA0);
    (*boundaryItr)->applyFieldListGhostBoundary(gradA);
    (*boundaryItr)->applyFieldListGhostBoundary(gradB);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Enforce boundary conditions on the fluid state Fields.
  // FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);

  FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CSPH, 0.0);
  FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CSPH, Vector::zero);
  FieldList<Dimension, Scalar> A0 = state.fields(HydroFieldNames::A0_CSPH, 0.0);
  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  FieldList<Dimension, Vector> gradA0 = state.fields(HydroFieldNames::gradA0_CSPH, Vector::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    // (*boundaryItr)->enforceFieldListBoundary(vol);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    if (compatibleEnergyEvolution()) (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
    (*boundaryItr)->enforceFieldListBoundary(m0);
    (*boundaryItr)->enforceFieldListBoundary(m1);
    (*boundaryItr)->enforceFieldListBoundary(A0);
    (*boundaryItr)->enforceFieldListBoundary(A);
    (*boundaryItr)->enforceFieldListBoundary(B);
    (*boundaryItr)->enforceFieldListBoundary(C);
    (*boundaryItr)->enforceFieldListBoundary(D);
    (*boundaryItr)->enforceFieldListBoundary(gradA0);
    (*boundaryItr)->enforceFieldListBoundary(gradA);
    (*boundaryItr)->enforceFieldListBoundary(gradB);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
dumpState(FileIO& file, string pathName) const {
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mVolume, pathName + "/volume");
  file.write(mM0, pathName + "/m0");
  file.write(mM1, pathName + "/m1");
  file.write(mM2, pathName + "/m2");
  file.write(mA0, pathName + "/A0");
  file.write(mA, pathName + "/A");
  file.write(mB, pathName + "/B");
  file.write(mC, pathName + "/C");
  file.write(mD, pathName + "/D");
  file.write(mGradA0, pathName + "/gradA0");
  file.write(mGradA, pathName + "/gradA");
  file.write(mGradB, pathName + "/gradB");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
restoreState(const FileIO& file, string pathName) {
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mVolume, pathName + "/volume");
  file.read(mM0, pathName + "/m0");
  file.read(mM1, pathName + "/m1");
  file.read(mM2, pathName + "/m2");
  file.read(mA0, pathName + "/A0");
  file.read(mA, pathName + "/A");
  file.read(mB, pathName + "/B");
  file.read(mC, pathName + "/C");
  file.read(mD, pathName + "/D");
  file.read(mGradA0, pathName + "/gradA0");
  file.read(mGradA, pathName + "/gradA");
  file.read(mGradB, pathName + "/gradB");
}

}
}

