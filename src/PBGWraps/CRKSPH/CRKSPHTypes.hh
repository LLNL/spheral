#ifndef __PBGWRAPS_CRKSPHTYPES__
#define __PBGWRAPS_CRKSPHTYPES__

#include "Geometry/Dimension.hh"
#include "CRKSPH/CRKSPHUtilities.hh"
#include "CRKSPH/CRKSPHHydroBase.hh"
#include "CRKSPH/CRKSPHHydroBaseRZ.hh"
#include "CRKSPH/SolidCRKSPHHydroBase.hh"
#include "CRKSPH/SolidCRKSPHHydroBaseRZ.hh"
#include "CRKSPH/CRKSPHVariant.hh"
#include "CRKSPH/computeVoronoiVolume.hh"
#include "CRKSPH/computeOccupancyVolume.hh"
#include "CRKSPH/computeCRKSPHSumVolume.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.hh"
#include "CRKSPH/computeSolidCRKSPHSumMassDensity.hh"
#include "CRKSPH/computeCRKSPHMoments.hh"
#include "CRKSPH/detectSurface.hh"
#include "CRKSPH/computeCRKSPHCorrections.hh"
#include "CRKSPH/centerOfMass.hh"
#include "CRKSPH/computeHullVolumes.hh"
#include "CRKSPH/computeNeighborHull.hh"
#include "CRKSPH/computeHVolumes.hh"
#include "CRKSPH/computeOccupancyVolume.hh"
#include "CRKSPH/gradientCRKSPH.hh"
#include "CRKSPH/interpolateCRKSPH.hh"
#include "SPH/NodeCoupling.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef CRKSPHHydroBase<Dim<1> > CRKSPHHydroBase1d;
typedef CRKSPHHydroBase<Dim<2> > CRKSPHHydroBase2d;
typedef CRKSPHHydroBase<Dim<3> > CRKSPHHydroBase3d;

typedef SolidCRKSPHHydroBase<Dim<1> > SolidCRKSPHHydroBase1d;
typedef SolidCRKSPHHydroBase<Dim<2> > SolidCRKSPHHydroBase2d;
typedef SolidCRKSPHHydroBase<Dim<3> > SolidCRKSPHHydroBase3d;

typedef CRKSPHVariant<Dim<1> > CRKSPHVariant1d;
typedef CRKSPHVariant<Dim<2> > CRKSPHVariant2d;
typedef CRKSPHVariant<Dim<3> > CRKSPHVariant3d;

//------------------------------------------------------------------------------
// Annoyingly we have to explicity disambiguate these.
//------------------------------------------------------------------------------
inline
void
CRKSPHKernelAndGradient1d(Dim<1>::Scalar* WCRKSPH,
                          Dim<1>::Scalar* gradWSPH,
                          Dim<1>::Vector& gradWCRKSPH,
                          const TableKernel<Dim<1> >& W,
                          const CRKOrder correctionOrder,
                          const Dim<1>::Vector& rij,
                          const Dim<1>::Vector& etai,
                          const Dim<1>::SymTensor& Hi,
                          const Dim<1>::Scalar& Hdeti,
                          const Dim<1>::Vector& etaj,
                          const Dim<1>::SymTensor& Hj,
                          const Dim<1>::Scalar& Hdetj,
                          const Dim<1>::Scalar& Ai,
                          const Dim<1>::Vector& Bi,
                          const Dim<1>::Tensor& Ci,
                          const Dim<1>::Vector& gradAi,
                          const Dim<1>::Tensor& gradBi,
                          const Dim<1>::ThirdRankTensor& gradCi,
                          const Dim<1>::Scalar correctionMin = std::numeric_limits<Dim<1>::Scalar>::lowest(),
                          const Dim<1>::Scalar correctionMax = std::numeric_limits<Dim<1>::Scalar>::max()) {
  return CRKSPHKernelAndGradient(*WCRKSPH, *gradWSPH, gradWCRKSPH, W, correctionOrder, rij, etai, Hi, Hdeti, etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, correctionMin, correctionMax);
}

inline
void
CRKSPHKernelAndGradient2d(Dim<2>::Scalar* WCRKSPH,
                          Dim<2>::Scalar* gradWSPH,
                          Dim<2>::Vector& gradWCRKSPH,
                          const TableKernel<Dim<2> >& W,
                          const CRKOrder correctionOrder,
                          const Dim<2>::Vector& rij,
                          const Dim<2>::Vector& etai,
                          const Dim<2>::SymTensor& Hi,
                          const Dim<2>::Scalar& Hdeti,
                          const Dim<2>::Vector& etaj,
                          const Dim<2>::SymTensor& Hj,
                          const Dim<2>::Scalar& Hdetj,
                          const Dim<2>::Scalar& Ai,
                          const Dim<2>::Vector& Bi,
                          const Dim<2>::Tensor& Ci,
                          const Dim<2>::Vector& gradAi,
                          const Dim<2>::Tensor& gradBi,
                          const Dim<2>::ThirdRankTensor& gradCi,
                          const Dim<2>::Scalar correctionMin = std::numeric_limits<Dim<2>::Scalar>::lowest(),
                          const Dim<2>::Scalar correctionMax = std::numeric_limits<Dim<2>::Scalar>::max()) {
  return CRKSPHKernelAndGradient(*WCRKSPH, *gradWSPH, gradWCRKSPH, W, correctionOrder, rij, etai, Hi, Hdeti, etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, correctionMin, correctionMax);
}

inline
void
CRKSPHKernelAndGradient3d(Dim<3>::Scalar* WCRKSPH,
                          Dim<3>::Scalar* gradWSPH,
                          Dim<3>::Vector& gradWCRKSPH,
                          const TableKernel<Dim<3> >& W,
                          const CRKOrder correctionOrder,
                          const Dim<3>::Vector& rij,
                          const Dim<3>::Vector& etai,
                          const Dim<3>::SymTensor& Hi,
                          const Dim<3>::Scalar& Hdeti,
                          const Dim<3>::Vector& etaj,
                          const Dim<3>::SymTensor& Hj,
                          const Dim<3>::Scalar& Hdetj,
                          const Dim<3>::Scalar& Ai,
                          const Dim<3>::Vector& Bi,
                          const Dim<3>::Tensor& Ci,
                          const Dim<3>::Vector& gradAi,
                          const Dim<3>::Tensor& gradBi,
                          const Dim<3>::ThirdRankTensor& gradCi,
                          const Dim<3>::Scalar correctionMin = std::numeric_limits<Dim<3>::Scalar>::lowest(),
                          const Dim<3>::Scalar correctionMax = std::numeric_limits<Dim<3>::Scalar>::max()) {
  return CRKSPHKernelAndGradient(*WCRKSPH, *gradWSPH, gradWCRKSPH, W, correctionOrder, rij, etai, Hi, Hdeti, etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, correctionMin, correctionMax);
}

// //------------------------------------------------------------------------------
// // compputeCRKSPHSumMassDensity with a std::vector<Boundary> rather than iterators.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// computeCRKSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
//                             const TableKernel<Dimension>& W,
//                             const FieldList<Dimension, typename Dimension::Vector>& position,
//                             const FieldList<Dimension, typename Dimension::Scalar>& mass,
//                             const FieldList<Dimension, typename Dimension::SymTensor>& H,
//                             const std::vector<Boundary<Dimension>*>& boundaries,
//                             FieldList<Dimension, typename Dimension::Scalar>& massDensity) {
//   computeCRKSPHSumMassDensity(connectivityMap, W, position, mass, H, 
//                             boundaries.begin(), boundaries.end(),
//                             massDensity);
// }

}

#endif
