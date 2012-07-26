//---------------------------------Spheral++----------------------------------//
// GenericHydro -- The base class for all Spheral++ hydro implementations.
//
// Created by JMO, Sat May 20 22:50:20 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_GenericHydro__
#define __Spheral_GenericHydro__

#include "Physics.hh"
#include "Geometry/Dimension.hh"

#ifndef __GCCXML__
#include "Kernel/TableKernel.hh"
#endif

namespace Spheral {
  namespace ArtificialViscositySpace {
    template<typename Dimension> class ArtificialViscosity;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
}

namespace Spheral {
namespace PhysicsSpace {

// Many hydro algorithms have these sorts of choices for the mass density and H.
enum MassDensityType {
  SumDensity = 0,
  RigorousSumDensity = 1,
  IntegrateDensity = 2,
  VoronoiCellDensity = 3,
  SumVoronoiCellDensity = 4,
};

enum HEvolutionType {
  IdealH = 0,
  IntegrateH = 1,
};

template<typename Dimension>
class GenericHydro: public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  GenericHydro(const KernelSpace::TableKernel<Dimension>& W,
               const KernelSpace::TableKernel<Dimension>& WPi,
               ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
               const double cfl,
               const bool useVelocityMagnitudeForDt);

  // Destructor.
  virtual ~GenericHydro();

  // We require all Physics packages to provide a method returning their vote
  // for the next time step.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Allow access to the artificial viscosity.
  ArtificialViscositySpace::ArtificialViscosity<Dimension>& artificialViscosity() const;

  // Access the stored interpolation kernels.
  const KernelSpace::TableKernel<Dimension>& kernel() const;
  const KernelSpace::TableKernel<Dimension>& PiKernel() const;

  // Also allow access to the CFL timestep safety criteria.
  Scalar cfl() const;
  void cfl(Scalar cfl);

  // Attribute to determine if the absolute magnitude of the velocity should
  // be used in determining the timestep.
  bool useVelocityMagnitudeForDt() const;
  void useVelocityMagnitudeForDt(bool x);

  // Return the cumulative neighboring statistics.
  int minMasterNeighbor() const;
  int maxMasterNeighbor() const;
  double averageMasterNeighbor() const;

  int minCoarseNeighbor() const;
  int maxCoarseNeighbor() const;
  double averageCoarseNeighbor() const;

  int minRefineNeighbor() const;
  int maxRefineNeighbor() const;
  double averageRefineNeighbor() const;

  int minActualNeighbor() const;
  int maxActualNeighbor() const;
  double averageActualNeighbor() const;

protected:
  //-------------------------- Protected Interface --------------------------//
  void updateMasterNeighborStats(int numMaster) const;
  void updateCoarseNeighborStats(int numNeighbor) const;
  void updateRefineNeighborStats(int numNeighbor) const;
  void updateActualNeighborStats(int numNeighbor) const;

private:
  //--------------------------- Private Interface ---------------------------//
  ArtificialViscositySpace::ArtificialViscosity<Dimension>& mArtificialViscosity;
  const KernelSpace::TableKernel<Dimension>& mKernel;
  const KernelSpace::TableKernel<Dimension>& mPiKernel;
  Scalar mCfl;
  bool mUseVelocityMagnitudeForDt;

  mutable int mMinMasterNeighbor, mMaxMasterNeighbor, mSumMasterNeighbor;
  mutable int mMinCoarseNeighbor, mMaxCoarseNeighbor, mSumCoarseNeighbor;
  mutable int mMinRefineNeighbor, mMaxRefineNeighbor, mSumRefineNeighbor;
  mutable int mMinActualNeighbor, mMaxActualNeighbor, mSumActualNeighbor;
  mutable int mNormMasterNeighbor;
  mutable int mNormCoarseNeighbor;
  mutable int mNormRefineNeighbor;
  mutable int mNormActualNeighbor;

  // Forbidden methods.
  GenericHydro();
  GenericHydro(const GenericHydro&);
  GenericHydro& operator=(const GenericHydro&);
};

}
}

#ifndef __GCCXML__
#include "GenericHydroInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class GenericHydro;
  }
}

#endif
