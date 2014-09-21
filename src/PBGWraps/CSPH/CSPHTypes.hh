#ifndef __PBGWRAPS_CSPHTYPES__
#define __PBGWRAPS_CSPHTYPES__

#include "Geometry/Dimension.hh"
#include "CSPH/CSPHUtilities.hh"
#include "CSPH/CSPHHydroBase.hh"
#include "CSPH/computeCSPHSumMassDensity.hh"
#include "CSPH/computeCSPHCorrections.hh"
#include "CSPH/centerOfMass.hh"
#include "CSPH/computeHullVolumes.hh"
#include "CSPH/computeHVolumes.hh"
#include "CSPH/gradientCSPH.hh"
#include "CSPH/interpolateCSPH.hh"

namespace Spheral {
namespace CSPHSpace {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef CSPHHydroBase<Dim<1> > CSPHHydroBase1d;
typedef CSPHHydroBase<Dim<2> > CSPHHydroBase2d;
typedef CSPHHydroBase<Dim<3> > CSPHHydroBase3d;

//------------------------------------------------------------------------------
// Annoyingly we have to explicity disambiguate these.
//------------------------------------------------------------------------------
inline
void
CSPHKernelAndGradient1d(const KernelSpace::TableKernel<Dim<1> >& W,
                        const Dim<1>::Vector& rij,
                        const Dim<1>::Vector& etaj,
                        const Dim<1>::SymTensor& Hj,
                        const Dim<1>::Scalar& Hdetj,
                        const Dim<1>::Scalar& Ai,
                        const Dim<1>::Vector& Bi,
                        const Dim<1>::Vector& gradAi,
                        const Dim<1>::Tensor& gradBi,
                        Dim<1>::Scalar* WCSPH,
                        Dim<1>::Scalar* gradWSPH,
                        Dim<1>::Vector& gradWCSPH) {
  return CSPHKernelAndGradient(W, rij, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, *WCSPH, *gradWSPH, gradWCSPH);
}

inline
void
CSPHKernelAndGradient2d(const KernelSpace::TableKernel<Dim<2> >& W,
                        const Dim<2>::Vector& rij,
                        const Dim<2>::Vector& etaj,
                        const Dim<2>::SymTensor& Hj,
                        const Dim<2>::Scalar& Hdetj,
                        const Dim<2>::Scalar& Ai,
                        const Dim<2>::Vector& Bi,
                        const Dim<2>::Vector& gradAi,
                        const Dim<2>::Tensor& gradBi,
                        Dim<2>::Scalar* WCSPH,
                        Dim<2>::Scalar* gradWSPH,
                        Dim<2>::Vector& gradWCSPH) {
  return CSPHKernelAndGradient(W, rij, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, *WCSPH, *gradWSPH, gradWCSPH);
}

inline
void
CSPHKernelAndGradient3d(const KernelSpace::TableKernel<Dim<3> >& W,
                        const Dim<3>::Vector& rij,
                        const Dim<3>::Vector& etaj,
                        const Dim<3>::SymTensor& Hj,
                        const Dim<3>::Scalar& Hdetj,
                        const Dim<3>::Scalar& Ai,
                        const Dim<3>::Vector& Bi,
                        const Dim<3>::Vector& gradAi,
                        const Dim<3>::Tensor& gradBi,
                        Dim<3>::Scalar* WCSPH,
                        Dim<3>::Scalar* gradWSPH,
                        Dim<3>::Vector& gradWCSPH) {
  return CSPHKernelAndGradient(W, rij, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, *WCSPH, *gradWSPH, gradWCSPH);
}

//------------------------------------------------------------------------------
// compputeCSPHSumMassDensity with a std::vector<Boundary> rather than iterators.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
computeCSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                          const KernelSpace::TableKernel<Dimension>& W,
                          const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                          const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                          const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& volume,
                          const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                          const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                          FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity) {
  computeCSPHSumMassDensity(connectivityMap, W, position, mass, volume, H, 
                            boundaries.begin(), boundaries.end(),
                            massDensity);
}

}
}

#endif
