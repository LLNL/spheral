//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBaseRZ -- The axisymmetric (RZ) SPH/ASPH solid material
//                        hydrodynamic package for Spheral++.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Mon May  9 11:01:51 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidSPHHydroBaseRZ_hh__
#define __Spheral_SolidSPHHydroBaseRZ_hh__

#include <float.h>
#include <string>

#include "SolidSPHHydroBase.hh"
#include "Geometry/Dimension.hh"

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

class SolidSPHHydroBaseRZ: public SolidSPHHydroBase<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  typedef Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidSPHHydroBaseRZ(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                      DataBase<Dimension>& dataBase,
                      ArtificialViscosity<Dimension>& Q,
                      const TableKernel<Dimension>& W,
                      const TableKernel<Dimension>& WPi,
                      const TableKernel<Dimension>& WGrad,
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
                      const bool strengthInDamage,
                      const Vector& xmin,
                      const Vector& xmax);

  // Destructor.
  virtual ~SolidSPHHydroBaseRZ();

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
  virtual std::string label() const override { return "SolidSPHHydroBaseRZ"; }

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  SolidSPHHydroBaseRZ();
  SolidSPHHydroBaseRZ(const SolidSPHHydroBaseRZ&);
  SolidSPHHydroBaseRZ& operator=(const SolidSPHHydroBaseRZ&);
};

}

#endif
