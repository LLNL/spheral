#ifndef __PBGWRAPS_CRKSPHTYPES__
#define __PBGWRAPS_CRKSPHTYPES__

#include "Geometry/Dimension.hh"
#include "CRKSPH/CRKSPHUtilities.hh"
#include "CRKSPH/CRKSPHHydroBase.hh"
#include "CRKSPH/SolidCRKSPHHydroBase.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.hh"
#include "CRKSPH/computeHullSumMassDensity.hh"
#include "CRKSPH/computeCRKSPHCorrections.hh"
#include "CRKSPH/centerOfMass.hh"
#include "CRKSPH/computeVoronoiCentroids.hh"
#include "CRKSPH/computeHullVolumes.hh"
#include "CRKSPH/computeNeighborHull.hh"
#include "CRKSPH/computeHVolumes.hh"
#include "CRKSPH/gradientCRKSPH.hh"
#include "CRKSPH/interpolateCRKSPH.hh"
#include "SolidSPH/NodeCoupling.hh"

namespace Spheral {
namespace CRKSPHSpace {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef CRKSPHHydroBase<Dim<1> > CRKSPHHydroBase1d;
typedef CRKSPHHydroBase<Dim<2> > CRKSPHHydroBase2d;
typedef CRKSPHHydroBase<Dim<3> > CRKSPHHydroBase3d;

typedef SolidCRKSPHHydroBase<Dim<1> > SolidCRKSPHHydroBase1d;
typedef SolidCRKSPHHydroBase<Dim<2> > SolidCRKSPHHydroBase2d;
typedef SolidCRKSPHHydroBase<Dim<3> > SolidCRKSPHHydroBase3d;

//------------------------------------------------------------------------------
// Annoyingly we have to explicity disambiguate these.
//------------------------------------------------------------------------------
inline
void
CRKSPHKernelAndGradient1d(const KernelSpace::TableKernel<Dim<1> >& W,
                        const Dim<1>::Vector& rij,
                        const Dim<1>::Vector& etai,
                        const Dim<1>::SymTensor& Hi,
                        const Dim<1>::Scalar& Hdeti,
                        const Dim<1>::Vector& etaj,
                        const Dim<1>::SymTensor& Hj,
                        const Dim<1>::Scalar& Hdetj,
                        const Dim<1>::Scalar& Ai,
                        const Dim<1>::Vector& Bi,
                        const Dim<1>::Vector& gradAi,
                        const Dim<1>::Tensor& gradBi,
                        Dim<1>::Scalar* WCRKSPH,
                        Dim<1>::Scalar* gradWSPH,
                        Dim<1>::Vector& gradWCRKSPH) {
  return CRKSPHKernelAndGradient(W, rij, etai, Hi, Hdeti, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, *WCRKSPH, *gradWSPH, gradWCRKSPH);
}

inline
void
CRKSPHKernelAndGradient2d(const KernelSpace::TableKernel<Dim<2> >& W,
                        const Dim<2>::Vector& rij,
                        const Dim<2>::Vector& etai,
                        const Dim<2>::SymTensor& Hi,
                        const Dim<2>::Scalar& Hdeti,
                        const Dim<2>::Vector& etaj,
                        const Dim<2>::SymTensor& Hj,
                        const Dim<2>::Scalar& Hdetj,
                        const Dim<2>::Scalar& Ai,
                        const Dim<2>::Vector& Bi,
                        const Dim<2>::Vector& gradAi,
                        const Dim<2>::Tensor& gradBi,
                        Dim<2>::Scalar* WCRKSPH,
                        Dim<2>::Scalar* gradWSPH,
                        Dim<2>::Vector& gradWCRKSPH) {
  return CRKSPHKernelAndGradient(W, rij, etai, Hi, Hdeti, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, *WCRKSPH, *gradWSPH, gradWCRKSPH);
}

inline
void
CRKSPHKernelAndGradient3d(const KernelSpace::TableKernel<Dim<3> >& W,
                        const Dim<3>::Vector& rij,
                        const Dim<3>::Vector& etai,
                        const Dim<3>::SymTensor& Hi,
                        const Dim<3>::Scalar& Hdeti,
                        const Dim<3>::Vector& etaj,
                        const Dim<3>::SymTensor& Hj,
                        const Dim<3>::Scalar& Hdetj,
                        const Dim<3>::Scalar& Ai,
                        const Dim<3>::Vector& Bi,
                        const Dim<3>::Vector& gradAi,
                        const Dim<3>::Tensor& gradBi,
                        Dim<3>::Scalar* WCRKSPH,
                        Dim<3>::Scalar* gradWSPH,
                        Dim<3>::Vector& gradWCRKSPH) {
  return CRKSPHKernelAndGradient(W, rij, etai, Hi, Hdeti, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, *WCRKSPH, *gradWSPH, gradWCRKSPH);
}

// //------------------------------------------------------------------------------
// // compputeCRKSPHSumMassDensity with a std::vector<Boundary> rather than iterators.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// computeCRKSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
//                           const KernelSpace::TableKernel<Dimension>& W,
//                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
//                           const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
//                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
//                           const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
//                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity) {
//   computeCRKSPHSumMassDensity(connectivityMap, W, position, mass, H, 
//                             boundaries.begin(), boundaries.end(),
//                             massDensity);
// }

}
}

#endif
