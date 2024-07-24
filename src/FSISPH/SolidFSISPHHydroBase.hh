//---------------------------------Spheral++----------------------------------//
// SolidFSISPHHydroBase -- SolidSPHHydro modified to better handle 
//                         multimaterial material problems with interfaces
//                         and large density ratios. 
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidFSISPHHydroBase_hh__
#define __Spheral_SolidFSISPHHydroBase_hh__

#include <float.h>
#include <string>
#include <vector>

#include "Physics/GenericHydro.hh"

namespace Spheral {

enum class InterfaceMethod {
  HLLCInterface = 0,
  ModulusInterface = 1,
  NoInterface = 2,
};

enum class KernelAveragingMethod {
  NeverAverageKernels = 0,
  AlwaysAverageKernels = 1,
  AverageInterfaceKernels = 2,
};

enum class FSIMassDensityMethod {
  FSISumMassDensity = 0,
  PressureCorrectSumMassDensity = 1,
  HWeightedSumMassDensity = 2,
};

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class SmoothingScaleBase;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class SlideSurface;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class SolidFSISPHHydroBase: public GenericHydro<Dimension> {

public:

  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidFSISPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                    DataBase<Dimension>& dataBase,
                    ArtificialViscosity<Dimension>& Q,
                    SlideSurface<Dimension>& slide,
                    const TableKernel<Dimension>& W,
                    const double cfl,
                    const double surfaceForceCoefficient,
                    const double densityStabilizationCoefficient,
                    const double specificThermalEnergyDiffusionCoefficient,
                    const double xsphCoefficient,
                    const InterfaceMethod interfaceMethod,
                    const KernelAveragingMethod kernelAveragingMethod,
                    const std::vector<int> sumDensityNodeLists,
                    const bool useVelocityMagnitudeForDt,
                    const bool compatibleEnergyEvolution,
                    const bool evolveTotalEnergy,
                    const bool linearCorrectGradients,
                    const bool planeStrain,
                    const double interfacePmin,
                    const double interfaceNeighborAngleThreshold,
                    const FSIMassDensityMethod densityUpdate,
                    const HEvolutionType HUpdate,
                    const double epsTensile,
                    const double nTensile,
                    const Vector& xmin,
                    const Vector& xmax);

  virtual ~SolidFSISPHHydroBase();

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  virtual
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  virtual 
  void preStepInitialize(const DataBase<Dimension>& dataBase, 
                               State<Dimension>& state,
                               StateDerivatives<Dimension>& derivs) override;
  virtual
  void initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivatives) const override;
  void firstDerivativesLoop(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivatives) const;
  void secondDerivativesLoop(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivatives) const;

  virtual 
  void finalizeDerivatives(const Scalar time, 
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase, 
                           const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) const override;

  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  void linearReconstruction(const typename Dimension::Vector& ri,
                            const typename Dimension::Vector& rj,
                            const typename Dimension::Scalar& yi,
                            const typename Dimension::Scalar& yj,
                            const typename Dimension::Vector& DyDxi,
                            const typename Dimension::Vector& DyDxj,
                                  typename Dimension::Scalar& ytildei,
                                  typename Dimension::Scalar& ytildej) const;


  const TableKernel<Dimension>& kernel() const;
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;
  SlideSurface<Dimension>& slideSurface() const;

  FSIMassDensityMethod densityUpdate() const;
  void densityUpdate(FSIMassDensityMethod type);

  HEvolutionType HEvolution() const;
  void HEvolution(HEvolutionType type);

  InterfaceMethod interfaceMethod() const;
  void interfaceMethod(InterfaceMethod method);

  KernelAveragingMethod kernelAveragingMethod() const;
  void kernelAveragingMethod(KernelAveragingMethod method);

  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(bool val);

  bool evolveTotalEnergy() const;
  void evolveTotalEnergy(bool val);

  bool linearCorrectGradients() const;
  void linearCorrectGradients(bool val);

  bool applySelectSumDensity() const;
  void applySelectSumDensity(bool x);

  bool planeStrain() const;
  void planeStrain(bool val);

  std::vector<int> sumDensityNodeLists() const;
  void sumDensityNodeLists(std::vector<int> x);

  double surfaceForceCoefficient() const;
  void surfaceForceCoefficient(double x);

  double densityStabilizationCoefficient() const;
  void densityStabilizationCoefficient(double x);

  double specificThermalEnergyDiffusionCoefficient() const;
  void specificThermalEnergyDiffusionCoefficient(double x);

  double interfacePmin() const;
  void interfacePmin(double val);

  double interfaceNeighborAngleThreshold() const;
  void interfaceNeighborAngleThreshold(double val);

  double xsphCoefficient() const;
  void xsphCoefficient(double x);

  Scalar epsilonTensile() const;
  void epsilonTensile(Scalar val);

  Scalar nTensile() const;
  void nTensile(Scalar val);

  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  const std::vector<Vector>&             pairAccelerations() const;
  const std::vector<Scalar>&             pairDepsDt() const;

  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    damagedPressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, Scalar>&    bulkModulus() const;
  const FieldList<Dimension, Scalar>&    shearModulus() const;
  const FieldList<Dimension, Scalar>&    yieldStrength() const;
  const FieldList<Dimension, Scalar>&    plasticStrain0() const;
  const FieldList<Dimension, Scalar>&    inverseEquivalentDeviatoricStress() const;
  const FieldList<Dimension, Scalar>&    volume() const;
  const FieldList<Dimension, Vector>&    DxDt() const;
  const FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldList<Dimension, Vector>&    DvDt() const;
  const FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const;
  const FieldList<Dimension, SymTensor>& DHDt() const;
  const FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldList<Dimension, Vector>&    DPDx() const;
  const FieldList<Dimension, Vector>&    DepsDx() const;
  const FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldList<Dimension, Tensor>&    internalDvDx() const;
  const FieldList<Dimension, Tensor>&    M() const;
  const FieldList<Dimension, Tensor>&    localM() const;
  const FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  const FieldList<Dimension, Scalar>&    effectiveViscousPressure() const;
  const FieldList<Dimension, Scalar>&    normalization() const;
  const FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldList<Dimension, SymTensor>& massSecondMoment() const;

  const FieldList<Dimension, int>& interfaceFlags() const;
  const FieldList<Dimension, Vector>& interfaceAreaVectors() const;
  const FieldList<Dimension, Vector>& interfaceNormals() const;
  const FieldList<Dimension, Scalar>& interfaceSmoothness() const;

  const FieldList<Dimension, int>& newInterfaceFlags() const;
  const FieldList<Dimension, Vector>& newInterfaceAreaVectors() const;
  const FieldList<Dimension, Vector>& newInterfaceNormals() const;
  const FieldList<Dimension, Scalar>& interfaceSmoothnessNormalization() const;
  const FieldList<Dimension, Scalar>& interfaceFraction() const;
  const FieldList<Dimension, Scalar>& newInterfaceSmoothness() const;
  const FieldList<Dimension, Scalar>& interfaceAngles() const;
  
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidFSISPHHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
 //****************************************************************************

private:
  const TableKernel<Dimension>& mKernel;
  const SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;
  SlideSurface<Dimension>& mSlideSurface;

  FSIMassDensityMethod mDensityUpdate;
  HEvolutionType mHEvolution;
  InterfaceMethod mInterfaceMethod;                   // switch for material interface method
  KernelAveragingMethod mKernelAveragingMethod;       // how do we handle our kernels?

  bool mCompatibleEnergyEvolution;
  bool mEvolveTotalEnergy; 
  bool mLinearCorrectGradients;
  bool mPlaneStrain;
  bool mApplySelectDensitySum;                        // switch for density sum
  std::vector<int> mSumDensityNodeLists;              // turn on density sum subset of nodeLists

  double mSurfaceForceCoefficient;                    // Monaghan 2013 force increase @ interface
  double mDensityStabilizationCoefficient;            // adjusts DvDx to stabilize rho
  double mSpecificThermalEnergyDiffusionCoefficient;  // controls diffusion of eps
  double mXSPHCoefficient;                            // controls amount of xsph-ing
  double mInterfacePmin;                              // min pressure across interfaces (similar to eos)
  double mInterfaceNeighborAngleThreshold;            // opening angle used to id surf/interface nodes

  Scalar mEpsTensile;
  Scalar mnTensile;

  Vector mxmin;
  Vector mxmax;

  std::vector<Vector> mPairAccelerations;              // store pairwise contribution to DvDt for compatible
  std::vector<Scalar> mPairDepsDt;                     // store pairwise contribution to DepsDt for compatible

  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mDamagedPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldList<Dimension, Scalar>    mBulkModulus;
  FieldList<Dimension, Scalar>    mShearModulus;
  FieldList<Dimension, Scalar>    mYieldStrength;
  FieldList<Dimension, Scalar>    mPlasticStrain0;
  FieldList<Dimension, Scalar>    mInverseEquivalentDeviatoricStress;
  FieldList<Dimension, Scalar>    mVolume;
  FieldList<Dimension, Vector>    mDxDt;
  FieldList<Dimension, Vector>    mXSPHDeltaV;
  FieldList<Dimension, Scalar>    mXSPHWeightSum;
  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldList<Dimension, SymTensor> mDHDt;
  FieldList<Dimension, SymTensor> mHideal;
  FieldList<Dimension, Vector>    mDPDx;
  FieldList<Dimension, Vector>    mDepsDx;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Tensor>    mInternalDvDx;
  FieldList<Dimension, Tensor>    mM;
  FieldList<Dimension, Tensor>    mLocalM;
  FieldList<Dimension, Scalar>    mMaxViscousPressure;
  FieldList<Dimension, Scalar>    mEffViscousPressure;
  FieldList<Dimension, Scalar>    mNormalization;
  FieldList<Dimension, Scalar>    mWeightedNeighborSum;
  FieldList<Dimension, SymTensor> mMassSecondMoment;

  FieldList<Dimension, int> mInterfaceFlags;                  // flags indicating interface type
  FieldList<Dimension, Vector> mInterfaceAreaVectors;         // interface area vectors that can be used for BCs
  FieldList<Dimension, Vector> mInterfaceNormals;             // surface normals between nodelists     
  FieldList<Dimension, Scalar> mInterfaceSmoothness;          // smoothness metric (0-1) 
  
  FieldList<Dimension, int> mNewInterfaceFlags;                   // flags indicating interface type
  FieldList<Dimension, Vector> mNewInterfaceAreaVectors;          // interface area vectors that can be used for BCs next time step
  FieldList<Dimension, Vector> mNewInterfaceNormals;              // surface normals between nodelists next time step    
  FieldList<Dimension, Scalar> mInterfaceSmoothnessNormalization; // normalization for a our smoothness metric
  FieldList<Dimension, Scalar> mInterfaceFraction;                // normalization for same material nodes
  FieldList<Dimension, Scalar> mNewInterfaceSmoothness;           // smoothness metric (0-1) next time step 
  FieldList<Dimension, Scalar> mInterfaceAngles;                  // check the angle for free-surface master nodes (type 2 -> type 3)

  // No default constructor, copying, or assignment.
  SolidFSISPHHydroBase();
  SolidFSISPHHydroBase(const SolidFSISPHHydroBase&);
  SolidFSISPHHydroBase& operator=(const SolidFSISPHHydroBase&);

protected:
  //--------------------------- Protected Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;
};

}


#include "SolidFSISPHHydroBaseInline.hh"

#endif
