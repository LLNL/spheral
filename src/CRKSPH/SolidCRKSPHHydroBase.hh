//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidCRKSPHHydroBase_hh__
#define __Spheral_SolidCRKSPHHydroBase_hh__

#include <float.h>
#include <string>

#include "CRKSPH/CRKSPHHydroBase.hh"

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

template<typename Dimension>
class SolidCRKSPHHydroBase: public CRKSPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidCRKSPHHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                       ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                       const KernelSpace::TableKernel<Dimension>& W,
                       const KernelSpace::TableKernel<Dimension>& WPi,
                       const double filter,
                       const double cfl,
                       const bool useVelocityMagnitudeForDt,
                       const bool compatibleEnergyEvolution,
                       const bool XSPH,
                       const PhysicsSpace::MassDensityType densityUpdate,
                       const PhysicsSpace::HEvolutionType HUpdate,
                       const CRKSPHSpace::CRKOrder correctionOrder,
                       const CRKSPHSpace::CRKVolumeType volumeType,
                       const double epsTensile,
                       const double nTensile);

  // Destructor.
  virtual ~SolidCRKSPHHydroBase();

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

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  void initialize(const Scalar time,
                  const Scalar dt,
                  const DataBaseSpace::DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
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
  const FieldSpace::FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>& bulkModulus() const;
  const FieldSpace::FieldList<Dimension, Scalar>& shearModulus() const;
  const FieldSpace::FieldList<Dimension, Scalar>& yieldStrength() const;
  const FieldSpace::FieldList<Dimension, Scalar>& plasticStrain0() const;
  const FieldSpace::FieldList<Dimension, int>& fragIDs() const;

  const FieldSpace::FieldList<Dimension, Scalar>&          Adamage() const;
  const FieldSpace::FieldList<Dimension, Vector>&          Bdamage() const;
  const FieldSpace::FieldList<Dimension, Tensor>&          Cdamage() const;
  const FieldSpace::FieldList<Dimension, Vector>&          gradAdamage() const;
  const FieldSpace::FieldList<Dimension, Tensor>&          gradBdamage() const;
  const FieldSpace::FieldList<Dimension, ThirdRankTensor>& gradCdamage() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SolidCRKSPHHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldSpace::FieldList<Dimension, Scalar> mBulkModulus;
  FieldSpace::FieldList<Dimension, Scalar> mShearModulus;
  FieldSpace::FieldList<Dimension, Scalar> mYieldStrength;
  FieldSpace::FieldList<Dimension, Scalar> mPlasticStrain0;
  FieldSpace::FieldList<Dimension, int> mFragIDs;

  FieldSpace::FieldList<Dimension, Scalar>          mAdamage;
  FieldSpace::FieldList<Dimension, Vector>          mBdamage;
  FieldSpace::FieldList<Dimension, Tensor>          mCdamage;
  FieldSpace::FieldList<Dimension, Vector>          mGradAdamage;
  FieldSpace::FieldList<Dimension, Tensor>          mGradBdamage;
  FieldSpace::FieldList<Dimension, ThirdRankTensor> mGradCdamage;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif

  // No default constructor, copying, or assignment.
  SolidCRKSPHHydroBase();
  SolidCRKSPHHydroBase(const SolidCRKSPHHydroBase&);
  SolidCRKSPHHydroBase& operator=(const SolidCRKSPHHydroBase&);
};

}
}

#include "SolidCRKSPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace SolidCRKSPHSpace {
    template<typename Dimension> class SolidCRKSPHHydroBase;
  }
}

#endif
