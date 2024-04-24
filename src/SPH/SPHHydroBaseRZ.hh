//---------------------------------Spheral++----------------------------------//
// SPHHydroBaseRZ -- An SPH/ASPH hydrodynamic package for Spheral++,
//                   specialized for 2D RZ (cylindrical) geometry.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Fri May  6 16:18:36 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPHHydroBaseRZ_hh__
#define __Spheral_SPHHydroBaseRZ_hh__

#include <string>

#include "SPHHydroBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SPHHydroBaseRZ: public SPHHydroBase<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  typedef Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SPHHydroBaseRZ(DataBase<Dimension>& dataBase,
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
                 const double epsTensile,
                 const double nTensile,
                 const Vector& xmin,
                 const Vector& xmax);

  // Destructor.
  virtual ~SPHHydroBaseRZ();

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
  virtual std::string label() const override { return "SPHHydroBaseRZ" ; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  SPHHydroBaseRZ();
  SPHHydroBaseRZ(const SPHHydroBaseRZ&);
  SPHHydroBaseRZ& operator=(const SPHHydroBaseRZ&);
};

}

#endif
