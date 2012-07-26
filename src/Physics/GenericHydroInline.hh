#include <limits.h>
#include <float.h>

#include "DBC.hh"

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// Access the artificial viscosity.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ArtificialViscositySpace::ArtificialViscosity<Dimension>&
GenericHydro<Dimension>::artificialViscosity() const {
  return mArtificialViscosity;
}

//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const KernelSpace::TableKernel<Dimension>&
GenericHydro<Dimension>::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// Access the kernel used for artificial viscosity gradients.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const KernelSpace::TableKernel<Dimension>&
GenericHydro<Dimension>::
PiKernel() const {
  return mPiKernel;
}

//------------------------------------------------------------------------------
// Access the CFL safety criteria.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GenericHydro<Dimension>::cfl() const {
  return mCfl;
}

template<typename Dimension>
inline
void
GenericHydro<Dimension>::
cfl(typename Dimension::Scalar cfl) {
  mCfl = cfl;
}

//------------------------------------------------------------------------------
// Flag to set whether or not to use the magnitude of the velocity to set the
// timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GenericHydro<Dimension>::useVelocityMagnitudeForDt() const {
  return mUseVelocityMagnitudeForDt;
}

template<typename Dimension>
inline
void
GenericHydro<Dimension>::
useVelocityMagnitudeForDt(bool x) {
  mUseVelocityMagnitudeForDt = x;
}

//------------------------------------------------------------------------------
// Return the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minMasterNeighbor() const {
  return mMinMasterNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxMasterNeighbor() const {
  return mMaxMasterNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageMasterNeighbor() const {
  return double(mSumMasterNeighbor)/(mNormMasterNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the coarse neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minCoarseNeighbor() const {
  return mMinCoarseNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxCoarseNeighbor() const {
  return mMaxCoarseNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageCoarseNeighbor() const {
  return double(mSumCoarseNeighbor)/(mNormCoarseNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the refine neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minRefineNeighbor() const {
  return mMinRefineNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxRefineNeighbor() const {
  return mMaxRefineNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageRefineNeighbor() const {
  return double(mSumRefineNeighbor)/(mNormRefineNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the actual neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GenericHydro<Dimension>::minActualNeighbor() const {
  return mMinActualNeighbor;
}

template<typename Dimension>
inline
int
GenericHydro<Dimension>::maxActualNeighbor() const {
  return mMaxActualNeighbor;
}

template<typename Dimension>
inline
double
GenericHydro<Dimension>::averageActualNeighbor() const {
  return double(mSumActualNeighbor)/(mNormActualNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Update the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateMasterNeighborStats(int numMaster) const {
  if (numMaster > 0) {
    mMinMasterNeighbor = std::min(mMinMasterNeighbor, numMaster);
    mMaxMasterNeighbor = std::max(mMaxMasterNeighbor, numMaster);
    mSumMasterNeighbor += numMaster;
    mNormMasterNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the coarse neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateCoarseNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinCoarseNeighbor = std::min(mMinCoarseNeighbor, numNeighbor);
    mMaxCoarseNeighbor = std::max(mMaxCoarseNeighbor, numNeighbor);
    mSumCoarseNeighbor += numNeighbor;
    mNormCoarseNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the refine neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateRefineNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinRefineNeighbor = std::min(mMinRefineNeighbor, numNeighbor);
    mMaxRefineNeighbor = std::max(mMaxRefineNeighbor, numNeighbor);
    mSumRefineNeighbor += numNeighbor;
    mNormRefineNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the actual neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GenericHydro<Dimension>::
updateActualNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinActualNeighbor = std::min(mMinActualNeighbor, numNeighbor);
    mMaxActualNeighbor = std::max(mMaxActualNeighbor, numNeighbor);
    mSumActualNeighbor += numNeighbor;
    mNormActualNeighbor++;
  }
}

}
}

