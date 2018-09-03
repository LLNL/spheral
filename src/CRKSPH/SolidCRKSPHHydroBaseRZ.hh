//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The CRKSPH/ACRKSPH solid material hydrodynamic
// package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Fri May 13 10:50:36 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidCRKSPHHydroBaseRZ_hh__
#define __Spheral_SolidCRKSPHHydroBaseRZ_hh__

#include "CRKSPH/SolidCRKSPHHydroBase.hh"

#include <float.h>
#include <string>

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

class SolidCRKSPHHydroBaseRZ: public SolidCRKSPHHydroBase<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Dimension::FourthRankTensor FourthRankTensor;
  typedef Dimension::FifthRankTensor FifthRankTensor;

  typedef Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidCRKSPHHydroBaseRZ(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                         ArtificialViscosity<Dimension>& Q,
                         const TableKernel<Dimension>& W,
                         const TableKernel<Dimension>& WPi,
                         const double filter,
                         const double cfl,
                         const bool useVelocityMagnitudeForDt,
                         const bool compatibleEnergyEvolution,
                         const bool evolveTotalEnergy,
                         const bool XSPH,
                         const MassDensityType densityUpdate,
                         const HEvolutionType HUpdate,
                         const CRKOrder correctionOrder,
                         const CRKVolumeType volumeType,
                         const double epsTensile,
                         const double nTensile,
                         const bool damageRelieveRubble);

  // Destructor.
  virtual ~SolidCRKSPHHydroBaseRZ();

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

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Finalize the hydro at the completion of an integration step.
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) override;
               
  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // The state field lists we're maintaining.
  // In the RZ case we have the (theta,theta) component of the deviatoric stress.
  const FieldList<Dimension, Scalar>& deviatoricStressTT() const;
  const FieldList<Dimension, Scalar>& DdeviatoricStressTTDt() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SolidCRKSPHHydroBaseRZ"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Some internal scratch fields.
  FieldList<Dimension, Scalar> mDeviatoricStressTT;
  FieldList<Dimension, Scalar> mDdeviatoricStressTTDt;

  // No default constructor, copying, or assignment.
  SolidCRKSPHHydroBaseRZ();
  SolidCRKSPHHydroBaseRZ(const SolidCRKSPHHydroBaseRZ&);
  SolidCRKSPHHydroBaseRZ& operator=(const SolidCRKSPHHydroBaseRZ&);
};

}

#include "SolidCRKSPHHydroBaseRZInline.hh"

#else

// Forward declaration.
namespace Spheral {
  class SolidCRKSPHHydroBaseRZ;
}

#endif
