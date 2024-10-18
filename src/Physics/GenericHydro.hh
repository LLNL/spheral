//---------------------------------Spheral++----------------------------------//
// GenericHydro -- The base class for all Spheral++ hydro implementations.
//
// Created by JMO, Sat May 20 22:50:20 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_GenericHydro__
#define __Spheral_GenericHydro__

#include "Physics.hh"
#include "Geometry/Dimension.hh"
#include "SmoothingScale/SmoothingScaleBase.hh"

namespace Spheral {

template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class DataBase;

// Many hydro algorithms have these sorts of choices for the mass density and H.
enum class MassDensityType {
  SumDensity = 0,
  RigorousSumDensity = 1,
  HybridSumDensity = 2,
  IntegrateDensity = 3,
  VoronoiCellDensity = 4,
  SumVoronoiCellDensity = 5,
  CorrectedSumDensity = 6,
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
  GenericHydro(ArtificialViscosity<Dimension>& Q,
               const double cfl,
               const bool useVelocityMagnitudeForDt);

  // Destructor.
  virtual ~GenericHydro();

  // We require all Physics packages to provide a method returning their vote
  // for the next time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Allow access to the artificial viscosity.
  ArtificialViscosity<Dimension>& artificialViscosity() const;

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

  // Stored attributes about the last timestep chosen
  size_t DTrank() const;
  size_t DTNodeList() const;
  size_t DTnode() const;
  std::string DTreason() const;

protected:
  //-------------------------- Protected Interface --------------------------//
  void updateMasterNeighborStats(int numMaster) const;
  void updateCoarseNeighborStats(int numNeighbor) const;
  void updateRefineNeighborStats(int numNeighbor) const;
  void updateActualNeighborStats(int numNeighbor) const;

private:
  //--------------------------- Private Interface ---------------------------//
  ArtificialViscosity<Dimension>& mArtificialViscosity;
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
  mutable size_t mDTrank, mDTNodeList, mDTnode;
  mutable std::string mDTreason;

  // Forbidden methods.
  GenericHydro();
  GenericHydro(const GenericHydro&);
  GenericHydro& operator=(const GenericHydro&);
};

}

#include "GenericHydroInline.hh"

#endif
