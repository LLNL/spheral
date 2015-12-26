//---------------------------------Spheral++----------------------------------//
// PSPHHydroBase -- The PSPH/APSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Wed Dec 16 20:52:02 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_PSPHHydroBase_hh__
#define __Spheral_PSPHHydroBase_hh__

#include <string>

#include "SPHHydroBase.hh"

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
namespace SPHSpace {

template<typename Dimension>
class PSPHHydroBase: public SPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  PSPHHydroBase(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
                const KernelSpace::TableKernel<Dimension>& W,
                const KernelSpace::TableKernel<Dimension>& WPi,
                const double filter,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool gradhCorrection,
                const bool XSPH,
                const bool HopkinsConductivity,
                const bool sumMassDensityOverAllNodeLists,
                const PhysicsSpace::MassDensityType densityUpdate,
                const PhysicsSpace::HEvolutionType HUpdate,
                const double epsTensile,
                const double nTensile,
                const Vector& xmin,
                const Vector& xmax);

  // Destructor.
  virtual ~PSPHHydroBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                     State<Dimension>& state);

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Post-state update. For PSPH this is where we recompute the PSPH pressure and corrections.
  virtual 
  void postStateUpdate(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivatives) const;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

  // Flag determining if we're applying Hopkins 2014 conductivity.
  bool HopkinsConductivity() const;
  void HopkinsConductivity(const bool val);

  // The state field lists we're maintaining.
  const FieldSpace::FieldList<Dimension, Scalar>&    gamma() const;
  const FieldSpace::FieldList<Dimension, Scalar>&    PSPHcorrection() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "PSPHHydroBase"; }
  virtual void dumpState(FileIOSpace::FileIO& file, std::string pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, std::string pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//
  bool mHopkinsConductivity;

  //PSPH Fields
  FieldSpace::FieldList<Dimension, Scalar>    mGamma;
  FieldSpace::FieldList<Dimension, Scalar>    mPSPHcorrection;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  PSPHHydroBase();
  PSPHHydroBase(const PSPHHydroBase&);
  PSPHHydroBase& operator=(const PSPHHydroBase&);
};

}
}

#include "PSPHHydroBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace PSPHSpace {
    template<typename Dimension> class PSPHHydroBase;
  }
}

#endif
