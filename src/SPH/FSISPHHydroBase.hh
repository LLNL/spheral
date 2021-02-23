//---------------------------------Spheral++----------------------------------//
// FSISPHHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Spheral_FSISPHHydroBase_hh__
#define __Spheral_FSISPHHydroBase_hh__

#include <float.h>
#include <string>
#include <vector>

#include "SPH/SPHHydroBase.hh"

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
class FSISPHHydroBase: public SPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  FSISPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
               DataBase<Dimension>& dataBase,
               ArtificialViscosity<Dimension>& Q,
               const TableKernel<Dimension>& W,
               const TableKernel<Dimension>& WPi,
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
               const Vector& xmin,
               const Vector& xmax);

  // Destructor.
  virtual ~FSISPHHydroBase();
                       
  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  double alpha() const;
  void alpha(double x);

  double diffusionCoefficient() const;
  void diffusionCoefficient(double x);

  std::vector<int> sumDensityNodeLists() const;
  void sumDensityNodeLists(std::vector<int> x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "FSISPHHydroBase"; }
 //****************************************************************************

private:
  double mAlpha;                               // generalized density exponent
  double mDiffusionCoefficient;                // controls diffusion of rho and eps
  std::vector<int> mSumDensityNodeLists;       // turn on density sum subset of nodeLists

  // No default constructor, copying, or assignment.
  FSISPHHydroBase();
  FSISPHHydroBase(const FSISPHHydroBase&);
  FSISPHHydroBase& operator=(const FSISPHHydroBase&);

};

}

#include "FSISPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class FSISPHHydroBase;
}

#endif
