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

#include <float.h>
#include <string>

#include "CRKSPH/SolidCRKSPHHydroBase.hh"

namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace NodeSpace {
    template<typename Dimension> class SmoothingScaleBase;
  }
  namespace ArtificialViscositySpace {
    template<typename Dimension> class ArtificialViscosity;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace CRKSPHSpace {

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

  typedef PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidCRKSPHHydroBaseRZ(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                         ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                         const KernelSpace::TableKernel<Dimension>& W,
                         const KernelSpace::TableKernel<Dimension>& WPi,
                         const double filter,
                         const double cfl,
                         const bool useVelocityMagnitudeForDt,
                         const bool compatibleEnergyEvolution,
                         const bool evolveTotalEnergy,
                         const bool XSPH,
                         const PhysicsSpace::MassDensityType densityUpdate,
                         const PhysicsSpace::HEvolutionType HUpdate,
                         const CRKSPHSpace::CRKOrder correctionOrder,
                         const CRKSPHSpace::CRKVolumeType volumeType,
                         const double epsTensile,
                         const double nTensile);

  // Destructor.
  virtual ~SolidCRKSPHHydroBaseRZ();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                     State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs);

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Finalize the hydro at the completion of an integration step.
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBaseSpace::DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs);
               
  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

  // The state field lists we're maintaining.
  // In the RZ case we have the (theta,theta) component of the deviatoric stress.
  const FieldSpace::FieldList<Dimension, Scalar>& deviatoricStressTT() const;
  const FieldSpace::FieldList<Dimension, Scalar>& DdeviatoricStressTTDt() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SolidCRKSPHHydroBaseRZ"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, Scalar> mDeviatoricStressTT;
  FieldSpace::FieldList<Dimension, Scalar> mDdeviatoricStressTTDt;

  // No default constructor, copying, or assignment.
  SolidCRKSPHHydroBaseRZ();
  SolidCRKSPHHydroBaseRZ(const SolidCRKSPHHydroBaseRZ&);
  SolidCRKSPHHydroBaseRZ& operator=(const SolidCRKSPHHydroBaseRZ&);
};

}
}

#include "SolidCRKSPHHydroBaseRZInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace SolidCRKSPHSpace {
    class SolidCRKSPHHydroBaseRZ;
  }
}

#endif
