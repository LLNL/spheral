//---------------------------------Spheral++----------------------------------//
// SphericalSPHAreaWeighted -- An SPH/ASPH hydrodynamic package for Spheral++,
//                             specialized for 1D Spherical (r) geometry.
//
// This implementation uses linear-weighting (the spherical version of
// area-weighting in RZ cylindrical coordinates).n
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Thu Jun 29 13:59:18 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalSPHAreaWeighted_hh__
#define __Spheral_SphericalSPHAreaWeighted_hh__

#include <string>

#include "SPHHydroBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SphericalSPHAreaWeighted: public SPHHydroBase<Dim<1>> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<1> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  typedef Physics<Dimension>::TimeStepType TimeStepType;
  typedef Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SphericalSPHAreaWeighted(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                           DataBase<Dimension>& dataBase,
                           ArtificialViscosity<Dimension>& Q,
                           const TableKernel<Dimension>& W,
                           const TableKernel<Dimension>& WPi,
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
  virtual ~SphericalSPHAreaWeighted();

  // Override the generic hydro timestep choice
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // This method is called once at the beginning of a timestep, after all state registration.
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

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;
               
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SphericalSPHAreaWeighted" ; }
  //****************************************************************************

  // We also have a funny self-Q term for interactions near the origin.
  double Qself() const;
  void Qself(const double x);

private:
  //--------------------------- Private Interface ---------------------------//
  double mQself;

  // No default constructor, copying, or assignment.
  SphericalSPHAreaWeighted();
  SphericalSPHAreaWeighted(const SphericalSPHAreaWeighted&);
  SphericalSPHAreaWeighted& operator=(const SphericalSPHAreaWeighted&);
};

}

#else

// Forward declaration.
namespace Spheral {
  class SphericalSPHAreaWeighted;
}

#endif
