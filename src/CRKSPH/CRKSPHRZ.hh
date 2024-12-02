//---------------------------------Spheral++----------------------------------//
// CRKSPHHydroBaseRZ -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Thu May 12 15:25:24 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPHHydroBaseRZ_hh__
#define __Spheral_CRKSPHHydroBaseRZ_hh__

#include <string>

#include "CRKSPHHydroBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;
}

namespace Spheral {

class CRKSPHHydroBaseRZ: public CRKSPHHydroBase<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using Tensor = Dimension::Tensor;
  using ThirdRankTensor = Dimension::ThirdRankTensor;
  using FourthRankTensor = Dimension::FourthRankTensor;
  using FifthRankTensor = Dimension::FifthRankTensor;
  using SymTensor = Dimension::SymTensor;
  using FacetedVolume = Dimension::FacetedVolume;

  using ConstBoundaryIterator = Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  CRKSPHHydroBaseRZ(DataBase<Dimension>& dataBase,
                    ArtificialViscosity<Dimension>& Q,
                    const RKOrder order,
                    const double filter,
                    const double cfl,
                    const bool useVelocityMagnitudeForDt,
                    const bool compatibleEnergyEvolution,
                    const bool evolveTotalEnergy,
                    const bool XSPH,
                    const MassDensityType densityUpdate,
                    const double epsTensile,
                    const double nTensile);

  // No default constructor, copying, or assignment.
  CRKSPHHydroBaseRZ() = delete;
  CRKSPHHydroBaseRZ(const CRKSPHHydroBaseRZ&) = delete;
  CRKSPHHydroBaseRZ& operator=(const CRKSPHHydroBaseRZ&) = delete;

  // Destructor.
  virtual ~CRKSPHHydroBaseRZ();

  // An optional hook to initialize once when the problem is starting up.
  // Typically this is used to size arrays once all the materials and NodeLists have
  // been created.  It is assumed after this method has been called it is safe to
  // call Physics::registerState for instance to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

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
  virtual std::string label() const override { return "CRKSPHHydroBaseRZ"; }
  //****************************************************************************
};

}

#endif
