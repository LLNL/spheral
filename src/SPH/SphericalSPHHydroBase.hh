//---------------------------------Spheral++----------------------------------//
// SphericalSPH -- An SPH/ASPH hydrodynamic package for Spheral++ specialized
//                 for 1D Spherical (r) geometry.
//
// Based on the algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Tue Dec 22 10:04:21 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalSPHHydroBase_hh__
#define __Spheral_SphericalSPHHydroBase_hh__

#include <utility>   // pair
#INCLUDE <string>

#include "SPHHydroBase.hh"
#include "Kernel/SphericalKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SphericalSPH: public SPHHydroBase<Dim<1>, SphericalSPH> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<1>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using Tensor = Dimension::Tensor;
  using SymTensor = Dimension::SymTensor;

  using PairAccelerationsType = PairwiseField<Dimension, std::pair<Vector, Vector>>;
  using ConstBoundaryIterator = Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SphericalSPHHydroBase(DataBase<Dimension>& dataBase,
                        ArtificialViscosity<Dimension>& Q,
                        const SphericalKernel& W,
                        const SphericalKernel& WPi,
                        const double cfl,
                        const bool useVelocityMagnitudeForDt,
                        const bool compatibleEnergyEvolution,
                        const bool evolveTotalEnergy,
                        const bool gradhCorrection,
                        const bool XSPH,
                        const bool correctVelocityGradient,
                        const bool sumMassDensityOverAllNodeLists,
                        const MassDensityType densityUpdate,
                        const double epsTensile,
                        const double nTensile,
                        const Vector& xmin,
                        const Vector& xmax);

  // No default constructor, copying, or assignment.
  SphericalSPHHydroBase() = delete;
  SphericalSPHHydroBase(const SphericalSPHHydroBase&) = delete;
  SphericalSPHHydroBase& operator=(const SphericalSPHHydroBase&) = delete;

  // Destructor.
  virtual ~SphericalSPHHydroBase();

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
  virtual std::string label() const override { return "SphericalSPHHydroBase" ; }
  //****************************************************************************

  // Access the stored interpolation kernels.
  // These hide the base class "kernel" methods which return vanilla TableKernels.
  const SphericalKernel& kernel() const;
  const SphericalKernel& PiKernel() const;

  // We also have a funny self-Q term for interactions near the origin.
  double Qself() const;
  void Qself(const double x);

private:
  //--------------------------- Private Interface ---------------------------//
  double mQself;

  // The specialized kernels
  const SphericalKernel& mKernel;
  const SphericalKernel& mPiKernel;
};

}

#endif
