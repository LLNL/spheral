//---------------------------------Spheral++----------------------------------//
// SolidFSISPHHydroBase -- 
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidFSISPHHydroBase_hh__
#define __Spheral_SolidFSISPHHydroBase_hh__

#include <float.h>
#include <string>
#include <vector>
#include "SPH/SolidSPHHydroBase.hh"

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
                    const TableKernel<Dimension>& W,
                    const TableKernel<Dimension>& WPi,
                    const TableKernel<Dimension>& WGrad,
                    const double alpha,
                    const double diffusionCoefficient,
                    const std::vector<int> sumDensityNodeLists,
                    const double filter,
                    const double cfl,
                    const bool useVelocityMagnitudeForDt,
                    const bool compatibleEnergyEvolution,
                    const bool evolveTotalEnergy,
                    const bool gradhCorrection,
                    const bool XSPH,
                    const bool correctVelocityGradient,
                    const bool sumMassDensityOverAllNodeLists,
                    const MassDensityType densityUpdate,
                    const HEvolutionType HUpdate,
                    const double epsTensile,
                    const double nTensile,
                    const bool damageRelieveRubble,
                    const bool negativePressureInDamage,
                    const bool strengthInDamage,
                    const Vector& xmin,
                    const Vector& xmax);

  // Destructor.
  virtual ~SolidFSISPHHydroBase();

  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  void computeFSISPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                                   const TableKernel<Dimension>& W,
                                   const FieldList<Dimension, typename Dimension::Vector>& position,
                                   const FieldList<Dimension, typename Dimension::Scalar>& mass,
                                   const FieldList<Dimension, typename Dimension::Scalar>& bulkModulus,
                                   const FieldList<Dimension, typename Dimension::Scalar>& pressure,
                                   const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                   FieldList<Dimension, typename Dimension::Scalar>& massDensity);
  
  double alpha() const;
  void alpha(double x);

  double diffusionCoefficient() const;
  void diffusionCoefficient(double x);

  std::vector<int> sumDensityNodeLists() const;
  void sumDensityNodeLists(std::vector<int> x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidFSISPHHydroBase"; }
 //****************************************************************************

private:
  double mAlpha;                               // generalized density exponent
  double mDiffusionCoefficient;                // controls diffusion of rho and eps
  std::vector<int> mSumDensityNodeLists;       // turn on density sum subset of nodeLists

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
