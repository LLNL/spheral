//---------------------------------Spheral++----------------------------------//
// A fake main script to construct the pieces of a Spheral simulation for 
// calling from assorted host codes that don't want to use python.
//
// An instance of this class essentially takes the place of writing Spheral
// python script.
//
// This is implemented as a singleton so the host can create one out of the
// ether and access it whenever needed.
// 
// Created by JMO, Thu Feb 28 2013
//----------------------------------------------------------------------------//
#include "boost/lexical_cast.hpp"

#include "SpheralPseudoScript2D.hh"
#include "SolidMaterial/LinearPolynomialEquationOfState.hh"
#include "SolidMaterial/NullStrength.hh"
#include "Kernel/BSplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "NodeList/SPHSmoothingScale.hh"
#include "NodeList/ASPHSmoothingScale.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"
#include "Distributed/VoronoiRedistributeNodes.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/iterateIdealH.hh"
#include "Utilities/globalNodeIDsInline.hh"
#include "Distributed/BoundingVolumeDistributedBoundary.hh"
#include "Distributed/NestedGridDistributedBoundary.hh"
#include "Boundary/ReflectingBoundary.hh"
#include "Field/Field.hh"

namespace Spheral {

using namespace std;

using Material::PhysicalConstants;
using SolidMaterial::SolidEquationOfState;
using SolidMaterial::StrengthModel;
using SolidMaterial::LinearPolynomialEquationOfState;
using SolidMaterial::NullStrength;
using SolidMaterial::SolidNodeList;
using KernelSpace::TableKernel;
using KernelSpace::BSplineKernel;
using KernelSpace::GaussianKernel;
using KernelSpace::PiGaussianKernel;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::SPHSmoothingScale;
using NodeSpace::ASPHSmoothingScale;
using NeighborSpace::TreeNeighbor;
using NeighborSpace::NestedGridNeighbor;
using NeighborSpace::ConnectivityMap;
using ArtificialViscositySpace::ArtificialViscosity;
using ArtificialViscositySpace::MonaghanGingoldViscosity;
using ArtificialViscositySpace::TensorMonaghanGingoldViscosity;
using SolidSPHSpace::SolidSPHHydroBase;
using DataBaseSpace::DataBase;
using IntegratorSpace::CheapSynchronousRK2;
using BoundarySpace::Boundary;
using BoundarySpace::ConstantBoundary;
using BoundarySpace::ReflectingBoundary;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Local utility functions.
//------------------------------------------------------------------------------
namespace {

//------------------------------------------------------------------------------
// Copy a set of C arrays to a FieldList.  Assume all elements are doubles, and
// the sizing between the array and the number of elements in the FieldList is
// correct.
//------------------------------------------------------------------------------
template<typename Dimension>
void
copyArrayToScalarFieldList(const double* array,
                           FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numElements();
    std::copy(&array[k], &array[k] + n, &(*fieldList[i]->begin()));
    k += n;
  }
}

//------------------------------------------------------------------------------
// Copy a set of C arrays to a FieldList.  Assume all elements are ints, and
// the sizing between the array and the number of elements in the FieldList is
// correct.
//------------------------------------------------------------------------------
template<typename Dimension>
void
copyArrayToIntFieldList(const int* array,
                        FieldList<Dimension, int>& fieldList) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numElements();
    std::copy(&array[k], &array[k] + n, &(*fieldList[i]->begin()));
    k += n;
  }
}

void
copyArrayToVectorFieldList(const double* x_array,
                           const double* y_array,
                           FieldList<Dim<2>, Dim<2>::Vector>& fieldList) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numElements();
    for (unsigned j = 0; j != n; ++j) {
      fieldList(i,j).x(x_array[k + j]);
      fieldList(i,j).y(y_array[k + j]);
    }
    k += n;
  }
}

void
copyArrayToTensorFieldList(const double* xx_array,
                           const double* xy_array,
                           const double* yx_array,
                           const double* yy_array,
                           FieldList<Dim<2>, Dim<2>::Tensor>& fieldList) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numElements();
    for (unsigned j = 0; j != n; ++j) {
      fieldList(i,j).xx(xx_array[k + j]);
      fieldList(i,j).xy(xy_array[k + j]);
      fieldList(i,j).yx(yx_array[k + j]);
      fieldList(i,j).yy(yy_array[k + j]);
    }
    k += n;
  }
}

void
copyArrayToSymTensorFieldList(const double* xx_array,
                              const double* xy_array,
                              const double* yy_array,
                              FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numElements();
    for (unsigned j = 0; j != n; ++j) {
      fieldList(i,j).xx(xx_array[k + j]);
      fieldList(i,j).xy(xy_array[k + j]);
      fieldList(i,j).yy(yy_array[k + j]);
    }
    k += n;
  }
}

//------------------------------------------------------------------------------
// Copy a FieldList to a set of C arrays.  In this case we only copy internal
// values back to the C arrays.
//------------------------------------------------------------------------------
template<typename Dimension>
void
copyScalarFieldListToArray(const FieldList<Dimension, typename Dimension::Scalar>& fieldList,
                           double* array) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numInternalElements();
    std::copy(fieldList[i]->begin(), fieldList[i]->begin() + n, &array[k]);
    k += n;
  }
}

//------------------------------------------------------------------------------
// Copy a FieldList to a set of C arrays.  In this case we only copy internal
// values back to the C arrays.
//------------------------------------------------------------------------------
template<typename Dimension>
void
copyIntFieldListToArray(const FieldList<Dimension, int>& fieldList,
                        int* array) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numInternalElements();
    std::copy(fieldList[i]->begin(), fieldList[i]->begin() + n, &array[k]);
    k += n;
  }
}

void
copyVectorFieldListToArray(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                           double* x_array,
                           double* y_array) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numInternalElements();
    for (unsigned j = 0; j != n; ++j) {
      x_array[k + j] = fieldList(i,j).x();
      y_array[k + j] = fieldList(i,j).y();
    }
    k += n;
  }
}

void
copyTensorFieldListToArray(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList,
                           double* xx_array,
                           double* xy_array,
                           double* yx_array,
                           double* yy_array) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numInternalElements();
    for (unsigned j = 0; j != n; ++j) {
      xx_array[k + j] = fieldList(i,j).xx();
      xy_array[k + j] = fieldList(i,j).xy();
      yx_array[k + j] = fieldList(i,j).yx();
      yy_array[k + j] = fieldList(i,j).yy();
    }
    k += n;
  }
}

void
copySymTensorFieldListToArray(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList,
                              double* xx_array,
                              double* xy_array,
                              double* yy_array) {
  const unsigned nfields = fieldList.numFields();
  unsigned k = 0;
  for (unsigned i = 0; i != nfields; ++i) {
    const unsigned n = fieldList[i]->numInternalElements();
    for (unsigned j = 0; j != n; ++j) {
      xx_array[k + j] = fieldList(i,j).xx();
      xy_array[k + j] = fieldList(i,j).xy();
      yy_array[k + j] = fieldList(i,j).yy();
    }
    k += n;
  }
}

}

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
SpheralPseudoScript2D&
SpheralPseudoScript2D::
instance() {
   if (mInstancePtr == 0) mInstancePtr = new SpheralPseudoScript2D;
   CHECK(mInstancePtr != 0);
   return *mInstancePtr;
}

//------------------------------------------------------------------------------
// initialize (called once at beginning of simulation).
//------------------------------------------------------------------------------
void
SpheralPseudoScript2D::
initialize(const bool ASPH,
           const bool XSPH,
           const bool compatibleEnergy,
           const bool vGradCorrection,
           const bool hGradCorrection,
           const bool sumMassDensity,
           const bool useVelocityDt,
           const bool ScalarQ,
           const bool addDistributedBoundary,
           const int kernelType,
           const int piKernelType,
           const int gradKernelType,
           const unsigned nmats,
           const double CFL,
           const double hmin,
           const double hmax,
           const double hminratio,
           const double nPerh,
           const double Clinear,
           const double Cquadratic,
           const double xmin_x,
           const double xmin_y,
           const double xmax_x,
           const double xmax_y) {

  // Build the min & max bounds.
  const Vector xmin(xmin_x, xmin_y), xmax(xmax_x, xmax_y);

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  // Create internal units (cm, gm, usec).
  me.mUnitsPtr = boost::shared_ptr<PhysicalConstants>(new PhysicalConstants(0.01,     // unit length (m)
                                                                            0.001,    // unit mass (kg)
                                                                            1.0e-6)); // unit time (sec)

  // Construct the stand-in fake EOS and strength model.  The host code will
  // actually fill in the state fields Spheral normally uses these for.
  me.mEOSptr = boost::shared_ptr<SolidEquationOfState<Dimension> >(new LinearPolynomialEquationOfState<Dimension>(1.0, 0.1, 10.0, 
                                                                                                                  1.0, 0.0, 0.0, 0.0,
                                                                                                                  1.0, 0.0, 0.0,
                                                                                                                  100.0, 
                                                                                                                  *me.mUnitsPtr,
                                                                                                                  0.0, 0.0, 1e100,
                                                                                                                  Material::ZeroPressure));
  me.mStrengthModelPtr = boost::shared_ptr<StrengthModel<Dimension> >(new NullStrength<Dimension>());

  // Build the interpolation kernels.
  if(kernelType == 0) {
    me.mKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(BSplineKernel<Dimension>(), 1000));
  }
  else if(kernelType == 1) {
    me.mKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(GaussianKernel<Dimension>(3.0), 1000));
  }
  else if(kernelType == 2) {
    me.mKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(PiGaussianKernel<Dimension>(7.0), 1000));
  }
  else {
    me.mKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(BSplineKernel<Dimension>(), 1000));
  }
  // Build the interpolation kernel for artificial viscosity.
  if(piKernelType == 0) {
    me.mPiKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(BSplineKernel<Dimension>(), 1000));
  }
  else if(piKernelType == 1) {
    me.mPiKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(GaussianKernel<Dimension>(3.0), 1000));
  }
  else if(piKernelType == 2) {
    me.mPiKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(PiGaussianKernel<Dimension>(7.0), 1000));
  }
  else {
    me.mPiKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(BSplineKernel<Dimension>(), 1000));
  }
  // Build the interpolation kernel for the velocity gradient.
  if(gradKernelType == 0) {
    me.mGradKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(BSplineKernel<Dimension>(), 1000));
  }
  else if(gradKernelType == 1) {
    me.mGradKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(GaussianKernel<Dimension>(3.0), 1000));
  }
  else if(gradKernelType == 2) {
    me.mGradKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(PiGaussianKernel<Dimension>(7.0), 1000));
  }
  else {
    me.mGradKernelPtr = boost::shared_ptr<TableKernel<Dimension> >(new TableKernel<Dimension>(BSplineKernel<Dimension>(), 1000));
  }

  // Construct the NodeLists for our materials.
  me.mNodeLists.reserve(nmats);
  me.mNeighbors.reserve(nmats);
  for (unsigned imat = 0; imat != nmats; ++imat) {
    const string name = "NodeList " + boost::lexical_cast<string>(imat);
    me.mNodeLists.push_back(new SolidNodeList<Dimension>(name,
                                                         *me.mEOSptr,
                                                         *me.mStrengthModelPtr,
                                                         0,
                                                         0,
                                                         hmin,
                                                         hmax,
                                                         hminratio,
                                                         nPerh,
                                                         500,      // maxNumNeighbors -- not currently used
                                                         0.0,      // rhoMin -- we depend on the host code
                                                         1e100));  // rhoMax -- we depend on the host code
    me.mNeighbors.push_back(new NestedGridNeighbor<Dimension>(me.mNodeLists.back(),
                                                              NeighborSpace::GatherScatter,
                                                              31,                            // numGridLevels
                                                              (xmax - xmin).maxElement(),    // topGridCellSize
                                                              Vector::zero,                  // origin
                                                              me.mKernelPtr->kernelExtent(),
                                                              1));                           // gridCellInfluenceRadius
    // me.mNeighbors.push_back(new TreeNeighbor<Dimension>(me.mNodeLists.back(),
    //                                                     NeighborSpace::GatherScatter,
    //                                                     me.mKernelPtr->kernelExtent(),
    //                                                     xmin,
    //                                                     xmax));
    me.mNodeLists.back().registerNeighbor(me.mNeighbors.back());
  }

  // Build the database and add our NodeLists.
  me.mDataBasePtr = boost::shared_ptr<DataBase<Dimension> >(new DataBase<Dimension>());
  for (unsigned imat = 0; imat != nmats; ++imat) {
    me.mDataBasePtr->appendNodeList(me.mNodeLists[imat]);
  }

  // Build the hydro physics objects.
  if (ASPH) {
    me.mSmoothingScaleMethodPtr = boost::shared_ptr<SmoothingScaleBase<Dimension> >(new ASPHSmoothingScale<Dimension>());
  } else {
    me.mSmoothingScaleMethodPtr = boost::shared_ptr<SmoothingScaleBase<Dimension> >(new SPHSmoothingScale<Dimension>());
  }
  if (ScalarQ) {
    me.mQptr = boost::shared_ptr<ArtificialViscosity<Dimension> >(new MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, false, false));
  } else {
    me.mQptr = boost::shared_ptr<ArtificialViscosity<Dimension> >(new TensorMonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic));
  }
  me.mQptr->epsilon2(0.01);
  me.mHydroPtr = boost::shared_ptr<SolidSPHHydroBase<Dimension> >(new SolidSPHHydroBase<Dimension>(*me.mSmoothingScaleMethodPtr,
                                                                                                   *me.mKernelPtr,
                                                                                                   *me.mPiKernelPtr,
                                                                                                   *me.mGradKernelPtr,
                                                                                                   *me.mQptr,
                                                                                                   0.0,                              // filter
                                                                                                   CFL,                              // cfl
                                                                                                   useVelocityDt,                    // useVelocityMagnitudeForDt
                                                                                                   compatibleEnergy,                 // compatibleEnergyEvolution
                                                                                                   hGradCorrection,                  // gradhCorrection
                                                                                                   XSPH,                             // XSPH
                                                                                                   vGradCorrection,                  // correctVelocityGradient
                                                                                                   sumMassDensity,                   // sumMassDensityOverAllNodeLists
                                                                                                   PhysicsSpace::RigorousSumDensity, // densityUpdate
                                                                                                   PhysicsSpace::IdealH,             // HUpdate
                                                                                                   0.0,                              // epsTensile
                                                                                                   4.0,                              // nTensile
                                                                                                   xmin,                             // xmin
                                                                                                   xmax));                           // xmax

  // Build a time integrator.  We're not going to use this to advance state,
  // but the other methods are useful.
  me.mIntegratorPtr = boost::shared_ptr<CheapSynchronousRK2<Dimension> >(new CheapSynchronousRK2<Dimension>(*me.mDataBasePtr));
  me.mIntegratorPtr->appendPhysicsPackage(*me.mHydroPtr);

  // Remember if we need to do parallel stuff ourself.
  me.mAddDistributedBoundary = addDistributedBoundary;

  // Do the one-time initialization work for our packages.
  me.mHydroPtr->initializeProblemStartup(*me.mDataBasePtr);

  // // If requested, add the appropriate DistributedBoundary for Spheral to handle
  // // distributed ghost nodes.
  // if (addDistributedBoundary) {
  //   me.mHydroPtr->appendBoundary(BoundarySpace::NestedGridDistributedBoundary<Dimension>::instance());
  // }
}

//------------------------------------------------------------------------------
// initializeStep -- to be called once at the beginning of a cycle.
// Returns:    the time step vote
// Arguments:  nintpermat : array indicating the number of internal nodes per material
//             npermat : array indicating the total number of nodes per material
//------------------------------------------------------------------------------
double
SpheralPseudoScript2D::
initializeStep(const unsigned* nintpermat,
               const unsigned* npermat,
               const double* mass,
               const double* massDensity,
               const double* position_x,
               const double* position_y,
               const double* specificThermalEnergy,
               const double* velocity_x,
               const double* velocity_y,
               const double* Hfield_xx,
               const double* Hfield_xy,
               const double* Hfield_yy,
               const double* pressure,
               const double* deviatoricStress_xx,
               const double* deviatoricStress_xy,
               const double* deviatoricStress_yy,
               const double* soundSpeed,
               const double* bulkModulus,
               const double* shearModulus,
               const double* yieldStrength,
               const double* plasticStrain) {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  // Check the input and set numbers of nodes.
  const unsigned nmats = me.mNodeLists.size();
  me.mNumInternalNodes = vector<unsigned>(nmats);
  me.mNumHostGhostNodes = vector<unsigned>(nmats);
  for (unsigned imat = 0; imat != nmats; ++imat) {
    VERIFY(nintpermat[imat] <= npermat[imat]);
    me.mNumInternalNodes[imat] = nintpermat[imat];
    me.mNumHostGhostNodes[imat] = npermat[imat] - nintpermat[imat];
    me.mNodeLists[imat].numInternalNodes(me.mNumInternalNodes[imat]);
    me.mNodeLists[imat].numGhostNodes(me.mNumHostGhostNodes[imat]);
  }

  // Prepare the state and such. Call registerState and registerDerivatives later!
  me.mStatePtr = boost::shared_ptr<State<Dimension> >(new State<Dimension>());
  me.mDerivsPtr = boost::shared_ptr<StateDerivatives<Dimension> >(new StateDerivatives<Dimension>());

  // Copy the given state into Spheral's structures.
  SpheralPseudoScript2D::updateState(mass, 
                                   massDensity,
                                   position_x, position_y, 
                                   specificThermalEnergy,
                                   velocity_x, velocity_y, 
                                   Hfield_xx, Hfield_xy, Hfield_yy, 
                                   pressure, 
                                   deviatoricStress_xx, deviatoricStress_xy, deviatoricStress_yy,
                                   soundSpeed, 
                                   bulkModulus, 
                                   shearModulus, 
                                   yieldStrength,
                                   plasticStrain);

  // Vote on a time step and return it.
  me.mIntegratorPtr->lastDt(1e10);
  const double dt = me.mIntegratorPtr->selectDt(1e-100, 1e100, *me.mStatePtr, *me.mDerivsPtr);
  return dt;
}

//------------------------------------------------------------------------------
// Just update the state without reinitializing everything for intermediate
// step estimates of the derivatives.
//------------------------------------------------------------------------------
void
SpheralPseudoScript2D::
updateState(const double* mass,
            const double* massDensity,
            const double* position_x,
            const double* position_y,
            const double* specificThermalEnergy,
            const double* velocity_x,
            const double* velocity_y,
            const double* Hfield_xx,
            const double* Hfield_xy,
            const double* Hfield_yy,
            const double* pressure,
            const double* deviatoricStress_xx,
            const double* deviatoricStress_xy,
            const double* deviatoricStress_yy,
            const double* soundSpeed,
            const double* bulkModulus,
            const double* shearModulus,
            const double* yieldStrength,
            const double* plasticStrain) {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();
  const unsigned nmats = me.mNodeLists.size();

  // Size the NodeLists.
  for (unsigned imat = 0; imat != nmats; ++imat) {
    me.mNodeLists[imat].numInternalNodes(me.mNumInternalNodes[imat]);
    me.mNodeLists[imat].numGhostNodes(me.mNumHostGhostNodes[imat]);
  }

  // register the Hydro fields.
  me.mHydroPtr->registerState(*me.mDataBasePtr,*me.mStatePtr);
  me.mHydroPtr->registerDerivatives(*me.mDataBasePtr,*me.mDerivsPtr);

  // Pull the state fields.
  FieldList<Dimension, Scalar> m = me.mStatePtr->fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> pos = me.mStatePtr->fields(HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> vel = me.mStatePtr->fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> rho = me.mStatePtr->fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> eps = me.mStatePtr->fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, SymTensor> H = me.mStatePtr->fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> P = me.mStatePtr->fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> cs = me.mStatePtr->fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, SymTensor> S = me.mStatePtr->fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> ps = me.mStatePtr->fields(SolidFieldNames::plasticStrain, 0.0);
  FieldList<Dimension, Scalar> K = me.mStatePtr->fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = me.mStatePtr->fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = me.mStatePtr->fields(SolidFieldNames::yieldStrength, 0.0);

  // Fill in the material properties.
  if(mass != NULL) {
    copyArrayToScalarFieldList(mass, m);
  }
  if(position_x != NULL && position_y != NULL) {
    copyArrayToVectorFieldList(position_x, position_y, pos);
  }
  if(velocity_x != NULL && velocity_y != NULL) {
    copyArrayToVectorFieldList(velocity_x, velocity_y, vel);
  }
  if(massDensity != NULL) {
    copyArrayToScalarFieldList(massDensity, rho);
  }
  if(specificThermalEnergy != NULL) {
    copyArrayToScalarFieldList(specificThermalEnergy, eps);
  }
  if(Hfield_xx != NULL && Hfield_xy != NULL && Hfield_yy != NULL) {
    copyArrayToSymTensorFieldList(Hfield_xx, Hfield_xy, Hfield_yy, H);
  }
  if(pressure != NULL) {
    copyArrayToScalarFieldList(pressure, P);
  }
  if(deviatoricStress_xx != NULL && deviatoricStress_xy != NULL && deviatoricStress_yy != NULL) {
    copyArrayToSymTensorFieldList(deviatoricStress_xx, deviatoricStress_xy, deviatoricStress_yy, S);
  }
  if(soundSpeed != NULL) {
    copyArrayToScalarFieldList(soundSpeed, cs);
  }
  if(bulkModulus != NULL) {
    copyArrayToScalarFieldList(bulkModulus, K);
  }
  if(shearModulus != NULL) {
    copyArrayToScalarFieldList(shearModulus, mu);
  }
  if(yieldStrength != NULL) {
    copyArrayToScalarFieldList(yieldStrength, Y);
  }
  if(plasticStrain != NULL) {
    copyArrayToScalarFieldList(plasticStrain, ps);
  }

  // Add host code boundaries
  me.mHydroPtr->clearBoundaries();
  for (unsigned i = 0; i<me.mHostCodeBoundaries.size(); ++i) {
    me.mHydroPtr->appendBoundary(*me.mHostCodeBoundaries[i]);
  }

  // If requested, add the appropriate DistributedBoundary for Spheral to handle
  // distributed ghost nodes.
  if (me.mAddDistributedBoundary) {
    me.mHydroPtr->appendBoundary(BoundarySpace::NestedGridDistributedBoundary<Dimension>::instance());
  }

  // Initialize the integrator.  This set's neighbor, connectivity, ghost nodes, etc.
  me.mIntegratorPtr->initialize(*me.mStatePtr, *me.mDerivsPtr);
}

//------------------------------------------------------------------------------
// evaluateDerivatives
//------------------------------------------------------------------------------
void
SpheralPseudoScript2D::
evaluateDerivatives(double* massDensitySum,
                    double* DmassDensityDt,
                    double* DvelocityDt_x,
                    double* DvelocityDt_y,
                    double* DxDt,
                    double* DyDt,
                    double* DspecificThermalEnergyDt,
                    double* DvelocityDx_xx,
                    double* DvelocityDx_xy,
                    double* DvelocityDx_yx,
                    double* DvelocityDx_yy,
                    double* DHfieldDt_xx,
                    double* DHfieldDt_xy,
                    double* DHfieldDt_yy,
                    double* HfieldIdeal_xx,
                    double* HfieldIdeal_xy,
                    double* HfieldIdeal_yy,
                    double* DdeviatoricStressDt_xx,
                    double* DdeviatoricStressDt_xy,
                    double* DdeviatoricStressDt_yy,
                    double* qpressure,
                    double* qwork) {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  // Zero out the stored derivatives.
  me.mDerivsPtr->Zero();

  // Have Spheral evaluate the fluid derivatives.
  me.mIntegratorPtr->evaluateDerivatives(0.0, 1.0, *me.mDataBasePtr, *me.mStatePtr, *me.mDerivsPtr);
  me.mIntegratorPtr->finalizeDerivatives(0.0, 1.0, *me.mDataBasePtr, *me.mStatePtr, *me.mDerivsPtr);
  me.mIntegratorPtr->finalize(0.0, 1.0, *me.mStatePtr, *me.mDerivsPtr);

  // Pull the individual derivative state fields.
  FieldList<Dimension, Scalar> rhoSum = me.mStatePtr->fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> DrhoDt = me.mDerivsPtr->fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = me.mDerivsPtr->fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Vector> DposDt = me.mDerivsPtr->fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = me.mDerivsPtr->fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = me.mDerivsPtr->fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, SymTensor> DHDt = me.mDerivsPtr->fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = me.mDerivsPtr->fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> DSDt = me.mDerivsPtr->fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> effViscousPressure = me.mDerivsPtr->fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  FieldList<Dimension, Scalar> viscousWork = me.mDerivsPtr->fields(HydroFieldNames::viscousWork, 0.0);

  // Copy the fluid derivative state to the arguments.
  copyScalarFieldListToArray(rhoSum, massDensitySum);
  copyScalarFieldListToArray(DrhoDt, DmassDensityDt);
  copyVectorFieldListToArray(DvDt, DvelocityDt_x, DvelocityDt_y);
  copyVectorFieldListToArray(DposDt, DxDt, DyDt);
  copyScalarFieldListToArray(DepsDt, DspecificThermalEnergyDt);
  copyTensorFieldListToArray(DvDx, DvelocityDx_xx, DvelocityDx_xy, DvelocityDx_yx, DvelocityDx_yy);
  copySymTensorFieldListToArray(DHDt, DHfieldDt_xx, DHfieldDt_xy, DHfieldDt_yy);
  copySymTensorFieldListToArray(Hideal, HfieldIdeal_xx, HfieldIdeal_xy, HfieldIdeal_yy);
  copySymTensorFieldListToArray(DSDt, DdeviatoricStressDt_xx, DdeviatoricStressDt_xy, 
                                                              DdeviatoricStressDt_yy);
  copyScalarFieldListToArray(effViscousPressure, qpressure);
  copyScalarFieldListToArray(viscousWork, qwork);
}

//------------------------------------------------------------------------------
// addBoundary
//------------------------------------------------------------------------------
void
SpheralPseudoScript2D::
addBoundary(char * bcname, double p0, double p1, double n0, double n1) {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();
  // Build the point and normal vectors
  const Vector point(p0,p1), normal(n0,n1);

  if(bcname == "symmetry" || bcname == "rigid") {
    me.mHostCodeBoundaries.push_back(boost::shared_ptr<ReflectingBoundary<Dimension> >(new ReflectingBoundary<Dimension>( (GeomPlane<Dimension>(point,normal)))));
  }

}

//------------------------------------------------------------------------------
// Provide a method of iterating the initial H tensors to something reasonable.
//------------------------------------------------------------------------------
void
SpheralPseudoScript2D::
iterateHfield(double* Hfield_xx,
              double* Hfield_xy,
              double* Hfield_yy,
              const int maxIterations,
              const double tolerance) {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  // Copy the input H info to the Spheral H field.
  FieldList<Dimension, SymTensor> H = me.mDataBasePtr->globalHfield();
  copyArrayToSymTensorFieldList(Hfield_xx, Hfield_xy, Hfield_yy, H);

  // We call the Utilities helper method to do this job.
  const vector<Boundary<Dimension>*>& bcs = me.mIntegratorPtr->uniqueBoundaryConditions();
  iterateIdealH(*me.mDataBasePtr, bcs, *me.mKernelPtr, *me.mSmoothingScaleMethodPtr, 
                maxIterations, tolerance);

  // Copy the internal H field back to the output.
  copySymTensorFieldListToArray(H, Hfield_xx, Hfield_xy, Hfield_yy);
}

//------------------------------------------------------------------------------
// Update the connectivity between nodes using Spheral's internal neighbor
// finding.
//------------------------------------------------------------------------------
void
SpheralPseudoScript2D::
updateConnectivity() {
  
  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  // We use the time integrators setGhostNodes method for this.  If we're not
  // using any of Spheral's boundary conditions this will only cause the 
  // connectivity to be recomputed.  However, if Spheral's boundary conditions
  // are being used this will also throw away and recreate the ghost nodes
  // from scratch.
  me.mIntegratorPtr->setGhostNodes();
}

//------------------------------------------------------------------------------
// Get the current connectivity.
//------------------------------------------------------------------------------
void
SpheralPseudoScript2D::
getConnectivity(int*** numNeighbors,
                int**** connectivity) {
  
  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();
  
  // Get the ConnectivityMap.
  const ConnectivityMap<Dimension>& cm = me.mDataBasePtr->connectivityMap();

  // Read the data out to the 4 deep C array (yipes!)
  const unsigned nmats = me.mNodeLists.size();
  CHECK(cm.nodeLists().size() == nmats);
  numNeighbors = new int**[nmats];
  connectivity = new int***[nmats];
  for (unsigned imat = 0; imat != nmats; ++imat) {
    const unsigned ni = me.mNodeLists[imat].numInternalNodes();
    numNeighbors[imat] = new int*[ni];
    connectivity[imat] = new int**[ni];
    for (unsigned i = 0; i != ni; ++i) {
      numNeighbors[imat][i] = new int[nmats];
      connectivity[imat][i] = new int*[nmats];
      const vector<vector<int> >& fullConnectivity = cm.connectivityForNode(imat, i);
      CHECK(fullConnectivity.size() == nmats);
      for (unsigned jmat = 0; jmat != nmats; ++jmat) {
        const unsigned nneigh = fullConnectivity[jmat].size();
        numNeighbors[imat][i][jmat] = nneigh;
        connectivity[imat][i][jmat] = new int[nneigh];
        std::copy(fullConnectivity[jmat].begin(), fullConnectivity[jmat].end(), connectivity[imat][i][jmat]);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Get the number of materials.
//------------------------------------------------------------------------------
int
SpheralPseudoScript2D::
getNumMaterials() {
  return SpheralPseudoScript2D::instance().mNodeLists.size();
}

//------------------------------------------------------------------------------
// Get the number of nodes per material.
//------------------------------------------------------------------------------
int*
SpheralPseudoScript2D::
getNumNodes() {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  const unsigned nmats = me.mNodeLists.size();
  int* result = new int[nmats];
  for (unsigned i = 0; i != nmats; ++i) {
    result[i] = me.mNodeLists[i].numNodes();
  }
}

//------------------------------------------------------------------------------
// Get the number of internal nodes per material.
//------------------------------------------------------------------------------
int*
SpheralPseudoScript2D::
getNumInternalNodes() {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  const unsigned nmats = me.mNodeLists.size();
  int* result = new int[nmats];
  for (unsigned i = 0; i != nmats; ++i) {
    result[i] = me.mNodeLists[i].numInternalNodes();
  }
}

//------------------------------------------------------------------------------
// Get the number of ghost nodes per material.
//------------------------------------------------------------------------------
int*
SpheralPseudoScript2D::
getNumGhostNodes() {

  // Get our instance.
  SpheralPseudoScript2D& me = SpheralPseudoScript2D::instance();

  const unsigned nmats = me.mNodeLists.size();
  int* result = new int[nmats];
  for (unsigned i = 0; i != nmats; ++i) {
    result[i] = me.mNodeLists[i].numGhostNodes();
  }
}

//------------------------------------------------------------------------------
// Constructor (private).
//------------------------------------------------------------------------------
SpheralPseudoScript2D::
SpheralPseudoScript2D() {
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
SpheralPseudoScript2D::
~SpheralPseudoScript2D() {
}

//------------------------------------------------------------------------------
// Initialize the static instance pointer.
//------------------------------------------------------------------------------
SpheralPseudoScript2D* SpheralPseudoScript2D::mInstancePtr = 0;

}
