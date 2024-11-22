//---------------------------------Spheral++----------------------------------//
// SolidSphericalSPHHydroBase -- The SPH/ASPH solid material SPH hydrodynamic
//                               specialized for 1D Spherical (r) geometry.
//
// Based on the algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Tue Apr 26 16:28:55 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidSphericalSPHHydroBase_hh__
#define __Spheral_SolidSphericalSPHHydroBase_hh__

#include <float.h>
#include <string>

#include "SPH/SolidSPHHydroBase.hh"
#include "Kernel/SphericalKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SolidSphericalSPHHydroBase: public SolidSPHHydroBase<Dim<1>> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<1>;
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SolidSphericalSPHHydroBase(DataBase<Dimension>& dataBase,
                             ArtificialViscosity<Dimension>& Q,
                             const SphericalKernel& W,
                             const SphericalKernel& WPi,
                             const SphericalKernel& WGrad,
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
                             const double epsTensile,
                             const double nTensile,
                             const bool damageRelieveRubble,
                             const bool strengthInDamage,
                             const Vector& xmin,
                             const Vector& xmax);

  // No default constructor, copying, or assignment.
  SolidSphericalSPHHydroBase() = delete;
  SolidSphericalSPHHydroBase(const SolidSphericalSPHHydroBase&) = delete;
  SolidSphericalSPHHydroBase& operator=(const SolidSphericalSPHHydroBase&) = delete;

  // Destructor.
  virtual ~SolidSphericalSPHHydroBase();

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

  // Access the stored interpolation kernels.
  // These hide the base class "kernel" methods which return vanilla TableKernels.
  const SphericalKernel& kernel() const;
  const SphericalKernel& PiKernel() const;
  const SphericalKernel& GradKernel() const;

  // We also have a funny self-Q term for interactions near the origin.
  double Qself() const;
  void Qself(const double x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidSphericalSPHHydroBase"; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  double mQself;

  // The specialized kernels
  const SphericalKernel& mKernel;
  const SphericalKernel& mPiKernel;
  const SphericalKernel& mGradKernel;
};

}

#endif
