#ifndef __PBGWRAPS_NODELISTTYPES__
#define __PBGWRAPS_NODELISTTYPES__

#include "Geometry/Dimension.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "NodeList/FixedSmoothingScale.hh"
#include "NodeList/SPHSmoothingScale.hh"
#include "NodeList/ASPHSmoothingScale.hh"
#include "NodeList/generateVoidNodes.hh"
#include "NodeList/nthNodalMoment.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {
namespace NodeSpace {

// //------------------------------------------------------------------------------
// // Provide a non-iterator based interface to generateVoidNodes.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// generateVoidNodes(const std::vector<NodeList<Dimension>*>& nodeLists,
//                   const MeshSpace::Mesh<Dimension>& mesh,
//                   const typename Dimension::Vector& xmin,
//                   const typename Dimension::Vector& xmax,
//                   NodeList<Dimension>& voidNodes,
//                   const double voidThreshold) {
//   generateVoidNodes(nodeLists.begin(), nodeLists.end(), mesh, xmin, xmax,
//                     voidNodes, voidThreshold);
// }

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to nthNodalMoment
//------------------------------------------------------------------------------
template<typename Dimension, unsigned moment>
inline
FieldSpace::FieldList<Dimension, typename MomentTraits<Dimension, moment>::Moment>
nthNodalMoment(const std::vector<NodeList<Dimension>*>& nodeLists,
               const KernelSpace::TableKernel<Dimension>& W,
               const bool renormalize) {
  return nthNodalMoment<Dimension, typename std::vector<NodeList<Dimension>*>::const_iterator, moment>
    (nodeLists.begin(), nodeLists.end(), W, renormalize);
}

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to zerothAndFirstNodalMoments
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
zerothAndFirstNodalMoments(const std::vector<NodeList<Dimension>*>& nodeLists,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const bool useGradientAsKernel,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& zerothMoment,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& firstMoment) {
  zerothAndFirstNodalMoments<Dimension, typename std::vector<NodeList<Dimension>*>::const_iterator>
    (nodeLists.begin(), nodeLists.end(), W, useGradientAsKernel, zerothMoment, firstMoment);
}

}
}

#endif
