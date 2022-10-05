//---------------------------------Spheral++----------------------------------//
// SolidFSISPHHydroBase -- 
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidFSISPHHydroBase_hh__
#define __Spheral_SolidFSISPHHydroBase_hh__

#include <float.h>
#include <string>
#include <vector>
#include "SPH/SolidSPHHydroBase.hh"

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
class SolidFSISPHHydroBase: public SolidSPHHydroBase<Dimension> {

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
                    const double filter,
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
                    const bool gradhCorrection,
                    const bool XSPH,
                    const bool correctVelocityGradient,
                    const MassDensityType densityUpdate,
                    const HEvolutionType HUpdate,
                    const double epsTensile,
                    const double nTensile,
                    const double interfacePmin,
                    const bool planeStrain,
                    const bool damageRelieveRubble,
                    const bool strengthInDamage,
                    const Vector& xmin,
                    const Vector& xmax);

  virtual ~SolidFSISPHHydroBase();

  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

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

  virtual 
  void finalizeDerivatives(const Scalar time, 
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase, 
                           const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) const override;

  void computeMCorrection(const typename Dimension::Scalar time,
                          const typename Dimension::Scalar dt,
                          const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                                StateDerivatives<Dimension>& derivatives) const;

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

  bool planeStrain() const;
  void planeStrain(bool x);

  double interfacePmin() const;
  void interfacePmin(double x);
  
  double surfaceForceCoefficient() const;
  void surfaceForceCoefficient(double x);

  double densityStabilizationCoefficient() const;
  void densityStabilizationCoefficient(double x);

  double specificThermalEnergyDiffusionCoefficient() const;
  void specificThermalEnergyDiffusionCoefficient(double x);

  double xsphCoefficient() const;
  void xsphCoefficient(double x);

  bool applySelectSumDensity() const;
  void applySelectSumDensity(bool x);

  std::vector<int> sumDensityNodeLists() const;
  void sumDensityNodeLists(std::vector<int> x);

  const std::vector<Scalar>& pairDepsDt() const;
  SlideSurface<Dimension>& slideSurface() const;

  InterfaceMethod interfaceMethod() const;
  void interfaceMethod(InterfaceMethod method);

  KernelAveragingMethod kernelAveragingMethod() const;
  void kernelAveragingMethod(KernelAveragingMethod method);

  const FieldList<Dimension, Scalar>& rawPressure() const;
  const FieldList<Dimension, Vector>& DPDx() const;
  const FieldList<Dimension, Vector>& DepsDx() const;
  const FieldList<Dimension, Vector>& interfaceNormals() const;
  const FieldList<Dimension, Scalar>& interfaceFraction() const;
  const FieldList<Dimension, Scalar>& interfaceSmoothness() const;
  const FieldList<Dimension, Vector>& newInterfaceNormals() const;
  const FieldList<Dimension, Vector>& smoothedInterfaceNormals() const;
  const FieldList<Dimension, Scalar>& newInterfaceFraction() const;
  const FieldList<Dimension, Scalar>& newInterfaceSmoothness() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidFSISPHHydroBase"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
 //****************************************************************************

private:
  bool mPlaneStrain;                                  // switch to update deviatoric stress according to plane-strain model
  SlideSurface<Dimension>& mSlideSurface;             // ref to the obj tracking slideSurfs between nodelists
  double mInterfacePmin;                              // minimum pressure allowed between material interfaces
  double mSurfaceForceCoefficient;                    // Monaghan 2013 force increase @ interface
  double mDensityStabilizationCoefficient;            // adjusts DvDx to stabilize rho
  double mSpecificThermalEnergyDiffusionCoefficient;  // controls diffusion of eps
  double mXSPHCoefficient;                            // controls amount of xsph-ing
  InterfaceMethod mInterfaceMethod;                   // switch for material interface method
  KernelAveragingMethod mKernelAveragingMethod;       // how do we handle our kernels?

  bool   mApplySelectDensitySum;                      // switch for density sum
  std::vector<int> mSumDensityNodeLists;              // turn on density sum subset of nodeLists
  
  std::vector<Scalar> mPairDepsDt;                     // store pairwise contribution to DepsDt for compatible
 
  FieldList<Dimension, Scalar> mRawPressure;                  // material interface normals
  FieldList<Dimension, Vector> mDPDx;                         // pressure gradient     
  FieldList<Dimension, Vector> mDepsDx;                       // specific thermal energy gradient    
  FieldList<Dimension, Vector> mInterfaceNormals;             // surface normals between nodelists     
  FieldList<Dimension, Scalar> mInterfaceFraction;            // fraction of dissimilar neighbor volume     
  FieldList<Dimension, Scalar> mInterfaceSmoothness;          // smoothness metric (0-1)    
  FieldList<Dimension, Vector> mNewInterfaceNormals;          // surface normals between nodelists next time step    
  FieldList<Dimension, Vector> mSmoothedInterfaceNormals;     // SPH interp of surface normal
  FieldList<Dimension, Scalar> mNewInterfaceFraction;         // fraction of dissimilar neighbor volume     
  FieldList<Dimension, Scalar> mNewInterfaceSmoothness;       // smoothness metric (0-1) next time step 

  // No default constructor, copying, or assignment.
  SolidFSISPHHydroBase();
  SolidFSISPHHydroBase(const SolidFSISPHHydroBase&);
  SolidFSISPHHydroBase& operator=(const SolidFSISPHHydroBase&);
};

}


#include "SolidFSISPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SolidFSISPHHydroBase;
}

#endif
