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
#ifndef __SpheralPseudoScript3D_hh__
#define __SpheralPseudoScript3D_hh__

#include "boost/shared_ptr.hpp"
#include "boost/ptr_container/ptr_vector.hpp"

#include "Material/PhysicalConstants.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Strength/SolidNodeList.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "SolidSPH/SolidSPHHydroBase.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Integrator/CheapSynchronousRK2.hh"
#include "Boundary/ConstantBoundary.hh"

namespace Spheral {

class SpheralPseudoScript3D {
public:
  //------------------------===== Public Interface =====-----------------------//
  typedef Dim<3> Dimension;
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  // Get the instance.
  static SpheralPseudoScript3D& instance();

  // initialize -- should be called once at the beginning of a simulation.
  static void initialize(const bool ASPH,
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
                         const double hminmaxratio,
                         const double nPerh,
                         const double Clinear,
                         const double Cquadratic,
                         const double xmin_x,
                         const double xmin_y,
                         const double xmin_z,
                         const double xmax_x,
                         const double xmax_y,
                         const double xmax_z);

  // initializeStep -- should be called once at the beginning of a cycle.
  static double initializeStep(const unsigned* nintpermat,
                               const unsigned* npermat,
                               const double* mass,
                               const double* massDensity,
                               const double* position_x,
                               const double* position_y,
                               const double* position_z,
                               const double* specificThermalEnergy,
                               const double* velocity_x,
                               const double* velocity_y,
                               const double* velocity_z,
                               const double* Hfield_xx,
                               const double* Hfield_xy,
                               const double* Hfield_xz,
                               const double* Hfield_yy,
                               const double* Hfield_yz,
                               const double* Hfield_zz,
                               const double* pressure,
                               const double* deviatoricStress_xx,
                               const double* deviatoricStress_xy,
                               const double* deviatoricStress_xz,
                               const double* deviatoricStress_yy,
                               const double* deviatoricStress_yz,
                               const double* deviatoricStress_zz,
                               const double* soundSpeed,
                               const double* bulkModulus,
                               const double* shearModulus,
                               const double* yieldStrength,
                               const double* plasticStrain);

  // updateState -- updates values of state fields without resizing.
  static void updateState(const double* mass,
                          const double* massDensity,
                          const double* position_x,
                          const double* position_y,
                          const double* position_z,
                          const double* specificThermalEnergy,
                          const double* velocity_x,
                          const double* velocity_y,
                          const double* velocity_z,
                          const double* Hfield_xx,
                          const double* Hfield_xy,
                          const double* Hfield_xz,
                          const double* Hfield_yy,
                          const double* Hfield_yz,
                          const double* Hfield_zz,
                          const double* pressure,
                          const double* deviatoricStress_xx,
                          const double* deviatoricStress_xy,
                          const double* deviatoricStress_xz,
                          const double* deviatoricStress_yy,
                          const double* deviatoricStress_yz,
                          const double* deviatoricStress_zz,
                          const double* soundSpeed,
                          const double* bulkModulus,
                          const double* shearModulus,
                          const double* yieldStrength,
                          const double* plasticStrain);

  // evaluateDerivatives -- computes the fluid time derivatives.
  static void evaluateDerivatives(double* massDensitySum,
                                  double* DmassDensityDt,
                                  double* DvelocityDt_x,
                                  double* DvelocityDt_y,
                                  double* DvelocityDt_z,
                                  double* DxDt,
                                  double* DyDt,
                                  double* DzDt,
                                  double* DspecificThermalEnergyDt,
                                  double* DvelocityDx_xx,
                                  double* DvelocityDx_xy,
                                  double* DvelocityDx_xz,
                                  double* DvelocityDx_yx,
                                  double* DvelocityDx_yy,
                                  double* DvelocityDx_yz,
                                  double* DvelocityDx_zx,
                                  double* DvelocityDx_zy,
                                  double* DvelocityDx_zz,
                                  double* DHfieldDt_xx,
                                  double* DHfieldDt_xy,
                                  double* DHfieldDt_xz,
                                  double* DHfieldDt_yy,
                                  double* DHfieldDt_yz,
                                  double* DHfieldDt_zz,
                                  double* HfieldIdeal_xx,
                                  double* HfieldIdeal_xy,
                                  double* HfieldIdeal_xz,
                                  double* HfieldIdeal_yy,
                                  double* HfieldIdeal_yz,
                                  double* HfieldIdeal_zz,
                                  double* DdeviatoricStressDt_xx,
                                  double* DdeviatoricStressDt_xy,
                                  double* DdeviatoricStressDt_xz,
                                  double* DdeviatoricStressDt_yy,
                                  double* DdeviatoricStressDt_yz,
                                  double* DdeviatoricStressDt_zz,
                                  double* qpressure,
                                  double* qwork);

  // Add boundary conditions
  static void addBoundary(char * bcname,
                          double p0, double p1, double p2,
                          double n0, double n1, double n2);

  // Take an initial guess for the H tensor field and optimize it.
  static void iterateHfield(double* Hfield_xx,
                            double* Hfield_xy,
                            double* Hfield_xz,
                            double* Hfield_yy,
                            double* Hfield_yz,
                            double* Hfield_zz,
                            const int maxIterations = 50,
                            const double tolerance = 1.0e-4);

  // Update the connectivity between nodes using Spheral's internal neighbor
  // finding.
  static void updateConnectivity();

  // Get the current connectivity information.
  static void getConnectivity(int*** numNeighbors,
                              int**** connectivity);

  // Get the number of materials.
  static int getNumMaterials();

  // Get the number of nodes in each material.
  static int* getNumNodes();
  static int* getNumInternalNodes();
  static int* getNumGhostNodes();

private:
  //------------------------===== Private Interface =====----------------------//
  // The one and only instance.
  static SpheralPseudoScript3D* mInstancePtr;

  // Numbers of nodes per material.
  std::vector<unsigned> mNumInternalNodes, mNumHostGhostNodes;

  // Flag as to whether we're doing the DistributedBoundary or not.
  bool mAddDistributedBoundary;

  // The material data.
  boost::shared_ptr<Material::PhysicalConstants> mUnitsPtr;
  boost::shared_ptr<SolidMaterial::SolidEquationOfState<Dimension> > mEOSptr;
  boost::shared_ptr<SolidMaterial::StrengthModel<Dimension> > mStrengthModelPtr;

  // The NodeList data.
  boost::ptr_vector<NeighborSpace::Neighbor<Dimension> > mNeighbors;
  boost::ptr_vector<SolidMaterial::SolidNodeList<Dimension> > mNodeLists;

  // Hydro bits.
  boost::shared_ptr<KernelSpace::TableKernel<Dimension> > mKernelPtr;
  boost::shared_ptr<KernelSpace::TableKernel<Dimension> > mPiKernelPtr;
  boost::shared_ptr<KernelSpace::TableKernel<Dimension> > mGradKernelPtr;
  boost::shared_ptr<NodeSpace::SmoothingScaleBase<Dimension> > mSmoothingScaleMethodPtr;
  boost::shared_ptr<ArtificialViscositySpace::ArtificialViscosity<Dimension> > mQptr;
  boost::shared_ptr<SolidSPHSpace::SolidSPHHydroBase<Dimension> > mHydroPtr;

  // Integrator and state.
  boost::shared_ptr<IntegratorSpace::CheapSynchronousRK2<Dimension> > mIntegratorPtr;
  boost::shared_ptr<DataBaseSpace::DataBase<Dimension> > mDataBasePtr;
  boost::shared_ptr<State<Dimension> > mStatePtr;
  boost::shared_ptr<StateDerivatives<Dimension> > mDerivsPtr;

  // A boundary to hold the host code values.
  std::vector<boost::shared_ptr<BoundarySpace::Boundary<Dimension> > > mHostCodeBoundaries;

  // No public constructors, destructor, or assignment.
  SpheralPseudoScript3D();
  SpheralPseudoScript3D(const SpheralPseudoScript3D&);
  SpheralPseudoScript3D& operator=(const SpheralPseudoScript3D&);
  ~SpheralPseudoScript3D();
};

}

#endif

