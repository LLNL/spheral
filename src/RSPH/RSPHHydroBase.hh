//---------------------------------Spheral++----------------------------------//
// RSPHHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Spheral_RSPHHydroBase_hh__
#define __Spheral_RSPHHydroBase_hh__

#include <float.h>
#include <string>
#include <vector>

#include "Physics/GenericHydro.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class SmoothingScaleBase;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class RSPHHydroBase: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  RSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
               DataBase<Dimension>& dataBase,
               ArtificialViscosity<Dimension>& Q,
               const TableKernel<Dimension>& W,
               const double cfl,
               const bool useVelocityMagnitudeForDt,
               const bool XSPH,
               const bool correctVelocityGradient,
               const HEvolutionType HUpdate,
               const double epsTensile,
               const double nTensile,
               const Vector& xmin,
               const Vector& xmax);

  // Destructor.
  virtual ~RSPHHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;   

  virtual void initialize(const Scalar time,
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

  void computeHLLCstate( const Vector& rij,
                         int nodeListi,
                         int nodeListj,
                         int i,
                         int j,
                         const Scalar& Pi,  
                         const Scalar& Pj,
                         const Scalar& rhoi, 
                         const Scalar& rhoj,
                         const Vector& vi,   
                         const Vector& vj,
                         const Scalar& ci,   
                         const Scalar& cj, 
                         Vector& vstar,
                         Scalar& Pstar) const;
  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;


  // Flag to select how we want to evolve the H tensor.
  // the continuity equation.
  HEvolutionType HEvolution() const;
  void HEvolution(HEvolutionType type);

  // Flag to determine if we're using the grad h correction.
  //bool gradhCorrection() const;
  //void gradhCorrection(bool val);

  // Flag to determine if we're using the XSPH algorithm.
  bool XSPH() const;
  void XSPH(bool val);

  // Flag to determine if we're applying the linear correction for the velocity gradient.
  bool correctVelocityGradient() const;
  void correctVelocityGradient(bool val);

  // Flag to determine if the sum density definition extends over neighbor NodeLists.
  //bool sumMassDensityOverAllNodeLists() const;
  //void sumMassDensityOverAllNodeLists(bool val);

  // Fraction of position filtering to apply.
  //double filter() const;
  //void filter(double val);

  // Parameters for the tensile correction force at small scales.
  Scalar epsilonTensile() const;
  void epsilonTensile(Scalar val);

  Scalar nTensile() const;
  void nTensile(Scalar val);

  // Optionally we can provide a bounding box for use generating the mesh
  // for the Voronoi mass density update.
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  // Access the stored interpolation kernels.
  const TableKernel<Dimension>& kernel() const;
  //const TableKernel<Dimension>& PiKernel() const;

  // The object defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod() const;

  // The state field lists we're maintaining.
  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  //const FieldList<Dimension, Scalar>&    volume() const;
  //const FieldList<Dimension, Scalar>&    omegaGradh() const;
  //const FieldList<Dimension, Scalar>&    specificThermalEnergy0() const;
  //const FieldList<Dimension, Scalar>&    entropy() const;
  const FieldList<Dimension, SymTensor>& Hideal() const;
  //const FieldList<Dimension, Scalar>&    maxViscousPressure() const;
  //const FieldList<Dimension, Scalar>&    effectiveViscousPressure() const;
  //const FieldList<Dimension, Scalar>&    massDensityCorrection() const;
  //const FieldList<Dimension, Scalar>&    viscousWork() const;
  //const FieldList<Dimension, Scalar>&    massDensitySum() const;
  //const FieldList<Dimension, Scalar>&    normalization() const;
  const FieldList<Dimension, Scalar>&    weightedNeighborSum() const;
  const FieldList<Dimension, SymTensor>& massSecondMoment() const;
  const FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldList<Dimension, Tensor>&    M() const;
  const FieldList<Dimension, Tensor>&    localM() const;
  const FieldList<Dimension, Vector>&    DxDt() const;
  const FieldList<Dimension, Vector>&    DvDt() const;
  const FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldList<Dimension, SymTensor>& DHDt() const;
  const FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldList<Dimension, Tensor>&    internalDvDx() const;
  const std::vector<Vector>&             pairAccelerations() const;

  double alpha() const;
  void alpha(double x);

  double diffusionCoefficient() const;
  void diffusionCoefficient(double x);

  std::vector<int> sumDensityNodeLists() const;
  void sumDensityNodeLists(std::vector<int> x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "RSPHHydroBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
 //****************************************************************************

protected:
// The interpolation kernels.
  const TableKernel<Dimension>& mKernel;

  // The method defining how we evolve smoothing scales.
  const SmoothingScaleBase<Dimension>& mSmoothingScaleMethod;

  // A bunch of switches.
  HEvolutionType mHEvolution;
  bool mXSPH, mCorrectVelocityGradient;

  // Tensile correction.
  Scalar mEpsTensile, mnTensile;

  // Optional bounding box for generating the mesh.
  Vector mxmin, mxmax;

  // Some internal scratch fields.
  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;
  //FieldList<Dimension, Scalar>    mOmegaGradh;
  //FieldList<Dimension, Scalar>    mSpecificThermalEnergy0;
  FieldList<Dimension, Scalar>    mEntropy;

  FieldList<Dimension, SymTensor> mHideal;
  //FieldList<Dimension, Scalar>    mMaxViscousPressure;
  //FieldList<Dimension, Scalar>    mEffViscousPressure;
  //FieldList<Dimension, Scalar>    mMassDensityCorrection;
  //FieldList<Dimension, Scalar>    mViscousWork;
  //FieldList<Dimension, Scalar>    mMassDensitySum;
  //FieldList<Dimension, Scalar>    mNormalization;

  FieldList<Dimension, Scalar>    mWeightedNeighborSum;
  FieldList<Dimension, SymTensor> mMassSecondMoment;

  FieldList<Dimension, Scalar>    mXSPHWeightSum;
  FieldList<Dimension, Vector>    mXSPHDeltaV;

  FieldList<Dimension, Vector>    mDxDt;
  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, SymTensor> mDHDt;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Tensor>    mInternalDvDx;
  FieldList<Dimension, Tensor>    mM;
  FieldList<Dimension, Tensor>    mLocalM;

  //FieldList<Dimension, Scalar>    mVolume;

  std::vector<Vector>             mPairAccelerations;

   // The restart registration.
  RestartRegistrationType mRestart;

  double mAlpha;                               // generalized density exponent
  double mDiffusionCoefficient;                // controls diffusion of rho and eps
  std::vector<int> mSumDensityNodeLists;       // turn on density sum subset of nodeLists

  //storing gradients for linear approx Reimann Problem
  FieldList<Dimension, Vector> mLastDrhoDx;
  FieldList<Dimension, Vector> mLastDcDx;
  FieldList<Dimension, Vector> mLastDpDx;
  FieldList<Dimension, Tensor> mLastDvDx;

private:
  // No default constructor, copying, or assignment.
  RSPHHydroBase();
  RSPHHydroBase(const RSPHHydroBase&);
  RSPHHydroBase& operator=(const RSPHHydroBase&);

};

}

#include "RSPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class RSPHHydroBase;
}

#endif
