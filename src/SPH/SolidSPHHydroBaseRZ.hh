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

class SolidSPHHydroBaseRZ: public SolidSPHHydroBase<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  typedef PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SolidSPHHydroBaseRZ(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
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
                      const bool damageRelieveRubble,
                      const Vector& xmin,
                      const Vector& xmax);

  // Destructor.
  virtual ~SolidSPHHydroBaseRZ();

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
  virtual std::string label() const { return "SolidSPHHydroBaseRZ"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, Scalar> mDeviatoricStressTT;
  FieldSpace::FieldList<Dimension, Scalar> mDdeviatoricStressTTDt;

  // No default constructor, copying, or assignment.
  SolidSPHHydroBaseRZ();
  SolidSPHHydroBaseRZ(const SolidSPHHydroBaseRZ&);
  SolidSPHHydroBaseRZ& operator=(const SolidSPHHydroBaseRZ&);
};

}
}

#ifndef __GCCXML__
#include "SolidSPHHydroBaseRZInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SPHSpace {
    class SolidSPHHydroBaseRZ;
  }
}

#endif
