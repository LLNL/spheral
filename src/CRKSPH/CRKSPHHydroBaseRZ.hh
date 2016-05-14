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

class CRKSPHHydroBaseRZ: public CRKSPHHydroBase<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Dimension::FourthRankTensor FourthRankTensor;
  typedef Dimension::FifthRankTensor FifthRankTensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::FacetedVolume FacetedVolume;

  typedef PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CRKSPHHydroBaseRZ(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
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
  virtual ~CRKSPHHydroBaseRZ();

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

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "CRKSPHHydroBaseRZ"; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  CRKSPHHydroBaseRZ();
  CRKSPHHydroBaseRZ(const CRKSPHHydroBaseRZ&);
  CRKSPHHydroBaseRZ& operator=(const CRKSPHHydroBaseRZ&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace CRKSPHSpace {
    class CRKSPHHydroBaseRZ;
  }
}

#endif
