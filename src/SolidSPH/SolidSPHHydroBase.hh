//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidSPHHydroBase_hh__
#define __Spheral_SolidSPHHydroBase_hh__

#include <float.h>
#include <string>

#include "SPH/SPHHydroBase.hh"

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
namespace SolidSPHSpace {

template<typename Dimension>
class SolidSPHHydroBase: public SPHSpace::SPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidSPHHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                    ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const KernelSpace::TableKernel<Dimension>& WPi,
                    const KernelSpace::TableKernel<Dimension>& WGrad,
                    const double filter,
                    const double cfl,
                    const bool useVelocityMagnitudeForDt,
                    const bool compatibleEnergyEvolution,
                    const bool evolveTotalEnergy,
                    const bool gradhCorrection,
                    const bool XSPH,
                    const bool correctVelocityGradient,
                    const bool sumMassDensityOverAllNodeLists,
                    const PhysicsSpace::MassDensityType densityUpdate,
                    const PhysicsSpace::HEvolutionType HUpdate,
                    const double epsTensile,
                    const double nTensile,
                    const Vector& xmin,
                    const Vector& xmax);

  // Destructor.
  virtual ~SolidSPHHydroBase();

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

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

  // Gradient kernel
  const KernelSpace::TableKernel<Dimension>& GradKernel() const;

  // The state field lists we're maintaining.
  const FieldSpace::FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>& bulkModulus() const;
  const FieldSpace::FieldList<Dimension, Scalar>& shearModulus() const;
  const FieldSpace::FieldList<Dimension, Scalar>& yieldStrength() const;
  const FieldSpace::FieldList<Dimension, Scalar>& plasticStrain0() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SolidSPHHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  // Gradient kernel
  const KernelSpace::TableKernel<Dimension>& mGradKernel;
  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldSpace::FieldList<Dimension, Scalar> mBulkModulus;
  FieldSpace::FieldList<Dimension, Scalar> mShearModulus;
  FieldSpace::FieldList<Dimension, Scalar> mYieldStrength;
  FieldSpace::FieldList<Dimension, Scalar> mPlasticStrain0;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif

  // No default constructor, copying, or assignment.
  SolidSPHHydroBase();
  SolidSPHHydroBase(const SolidSPHHydroBase&);
  SolidSPHHydroBase& operator=(const SolidSPHHydroBase&);
};

}
}

#ifndef __GCCXML__
#include "SolidSPHHydroBaseInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SolidSPHSpace {
    template<typename Dimension> class SolidSPHHydroBase;
  }
}

#endif
