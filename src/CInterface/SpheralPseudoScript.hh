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
// This is the base class for SpheralPsuedoScript2D and SpheralPsuedoScript3D,
// the idea being to put the common code here to be shared by the two
// concrete dimensional types.
// 
// Created by JMO, Thu Feb 28 2013
//----------------------------------------------------------------------------//
#ifndef __SpheralPseudoScript_hh__
#define __SpheralPseudoScript_hh__

#include "Material/PhysicalConstants.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "NodeList/SolidNodeList.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Physics/Physics.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Integrator/CheapSynchronousRK2.hh"

namespace Spheral {

template<typename Dimension>
class SpheralPseudoScript {
public:
  //------------------------===== Public Interface =====-----------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef Field<Dimension, std::vector<double>> FlawStorageType;

  // Get the instance.
  static SpheralPseudoScript& instance();

  // initialize -- should be called once at the beginning of a simulation.
  static void initialize(const bool     RZ,
                         const bool     CRK,
                         const bool     ASPH,
                         const bool     XSPH,
                         const bool     compatibleEnergy,
                         const bool     totalEnergy,
                         const bool     vGradCorrection,
                         const bool     hGradCorrection,
                         const bool     sumMassDensity,
                         const bool     useVelocityDt,
                         const bool     ScalarQ,
                         const int      distributedBoundary,
                         const int      kernelType,
                         const int      piKernelType,
                         const int      gradKernelType,
                         const int      nbspline,
                         const int      crkorder,
                         const int      damage,
                         const unsigned nmats,
                         const double   CFL,
                         const double   hmin,
                         const double   hmax,
                         const double   hminmaxratio,
                         const double   nPerh,
                         const double   Clinear,
                         const double   Cquadratic,
                         const Vector&  xmin,
                         const Vector&  xmax);

  // initializeStep -- should be called once at the beginning of a cycle.
  static double initializeStep(const unsigned* nintpermat,
                               const unsigned* npermat,
                               const double*   mass,
                               const double*   massDensity,
                               const double**  position,
                               const double*   specificThermalEnergy,
                               const double**  velocity,
                               const double**  Hfield,
                               const double*   pressure,
                               const double**  deviatoricStress,
                               const double*   deviatoricStressTT,
                               const double*   soundSpeed,
                               const double*   bulkModulus,
                               const double*   shearModulus,
                               const double*   yieldStrength,
                               const double*   plasticStrain,
                               const double*   scalarDamage,
                               const int*      particleType);

  // updateState -- updates values of state fields without resizing.
  static void updateState(const double*  mass,
                          const double*  massDensity,
                          const double** positionx,
                          const double*  specificThermalEnergy,
                          const double** velocity,
                          const double** Hfield,
                          const double*  pressure,
                          const double** deviatoricStress,
                          const double*  deviatoricStressTT,
                          const double*  soundSpeed,
                          const double*  bulkModulus,
                          const double*  shearModulus,
                          const double*  yieldStrength,
                          const double*  plasticStrain,
                          const double*  scalarDamage,
                          const int*     particleType);

  // evaluateDerivatives -- computes the fluid time derivatives.
  static void evaluateDerivatives(double*  massDensitySum,
                                  double*  DmassDensityDt,
                                  double** DvelocityDt,
                                  double** DxDt,
                                  double*  DspecificThermalEnergyDt,
                                  double** DvelocityDx,
                                  double** DHfieldDt,
                                  double** HfieldIdeal,
                                  double** DdeviatoricStressDt,
                                  double*  DdeviatoricStressDtTT,
                                  double*  qpressure,
                                  double*  qwork);

  // Add boundary conditions
  static void addBoundary(const Vector& point,
                          const Vector& normal);

  // Add boundary conditions
  static void addPeriodicBoundary(const Vector& point1,
                                  const Vector& normal1,
                                  const Vector& point2,
                                  const Vector& normal2);

  // Take an initial guess for the H tensor field and optimize it.
  static void iterateHfield(double**     Hfield,
                            const int    maxIterations = 50,
                            const double tolerance = 1.0e-4);

  // Compute the fragment identification field.
  static void computeFragmentID(double* damage,
                                double  fragRadius,
                                double  fragDensity,
                                double  fragDamage,
                                int*    fragments);

  // Sample the SPH state to a structured mesh.
  static void sampleLatticeMesh(const Vector&  xmin,
                                const Vector&  xmax,
                                const int*     nsamples,
                                double*        latticeDensity);

  static void polyhedralMesh(int*           nnodes,
                             int*           nfaces,
                             int*           ncells,
                             double**       coords,
                             int**          facetonodes,
                             int**          celltofaces);

  static void fillVolume(const int*     nnodes,
                         const int*     nfaces,
                         const double** coords,
                         const int*     conn,
                         const double   spacing,
                         const int      domain,
                         const int      ndomains,
                         double*        volume,
                         int*           nparticles,
                         double**       sphcoords);

  static void generateCylFromRZ(const int*     nnodes,
                                const double** coords,
                                const double** htensor,
                                const double** volume,
                                const double   frac,
                                int*           nparticles,
                                double**       sphcoords,
                                double**       sphhtensor,
                                double**       sphvolume);

  // Update the connectivity between nodes using Spheral's internal neighbor
  // finding.
  static void updateConnectivity();

  // Get the current connectivity information.
  static void getConnectivity(int***  numNeighbors,
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
  static SpheralPseudoScript* mInstancePtr;

  // Numbers of nodes per material.
  std::vector<unsigned> mNumInternalNodes, mNumHostGhostNodes;

  // Damage flag
  bool mDamage;

  // Flag as to whether we're doing the DistributedBoundary or not.
  int mDistributedBoundary;

  // The material data.
  std::shared_ptr<PhysicalConstants> mUnitsPtr;
  std::shared_ptr<SolidEquationOfState<Dimension> > mEOSptr;
  std::shared_ptr<StrengthModel<Dimension> > mStrengthModelPtr;

  // The NodeList data.
  std::vector<std::shared_ptr<Neighbor<Dimension> > > mNeighbors;
  std::vector<std::shared_ptr<SolidNodeList<Dimension> > > mNodeLists;

  // Hydro bits.
  std::shared_ptr<TableKernel<Dimension> > mKernelPtr;
  std::shared_ptr<TableKernel<Dimension> > mPiKernelPtr;
  std::shared_ptr<TableKernel<Dimension> > mGradKernelPtr;
  std::shared_ptr<SmoothingScaleBase<Dimension> > mSmoothingScaleMethodPtr;
  std::shared_ptr<ArtificialViscosity<Dimension> > mQptr;
  std::shared_ptr<Physics<Dimension> > mHydroPtr;

  // Integrator and state.
  std::shared_ptr<CheapSynchronousRK2<Dimension> > mIntegratorPtr;
  std::shared_ptr<DataBase<Dimension> > mDataBasePtr;
  std::shared_ptr<State<Dimension>> mStatePtr;
  std::shared_ptr<StateDerivatives<Dimension>> mDerivsPtr;

  // A boundary to hold the host code values.
  std::vector<std::shared_ptr<Boundary<Dimension> > > mHostCodeBoundaries;

  // No public constructors, destructor, or assignment.
  SpheralPseudoScript();
  SpheralPseudoScript(const SpheralPseudoScript&);
  SpheralPseudoScript& operator=(const SpheralPseudoScript&);
  ~SpheralPseudoScript();
};

}

#endif

