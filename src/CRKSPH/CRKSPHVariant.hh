//---------------------------------Spheral++----------------------------------//
// CRKSPHVariant -- A development variant of CRKSPH for experimentation.
//
// Created by JMO, Thu Oct 12 14:24:43 PDT 2017
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPHVariant_hh__
#define __Spheral_CRKSPHVariant_hh__

#include <string>

#include "Physics/GenericHydro.hh"
#include "CRKSPHHydroBase.hh"

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
class CRKSPHVariant: public CRKSPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  typedef typename PhysicsSpace::Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CRKSPHVariant(const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
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
                const bool detectSurfaces,
                const double detectThreshold,
                const double sweepAngle,
                const double detectRange,
                const double epsTensile,
                const double nTensile);

  // Destructor.
  virtual ~CRKSPHVariant() override;

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  CRKSPHVariant();
  CRKSPHVariant(const CRKSPHVariant&);
  CRKSPHVariant& operator=(const CRKSPHVariant&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace CRKSPHSpace {
    template<typename Dimension> class CRKSPHVariant;
  }
}

#endif
