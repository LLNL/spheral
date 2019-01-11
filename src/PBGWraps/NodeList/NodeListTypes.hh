#ifndef __PBGWRAPS_NODELISTTYPES__
#define __PBGWRAPS_NODELISTTYPES__

#include "Geometry/Dimension.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "NodeList/SolidNodeList.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "NodeList/FixedSmoothingScale.hh"
#include "NodeList/SPHSmoothingScale.hh"
#include "NodeList/ASPHSmoothingScale.hh"
#include "NodeList/generateVoidNodes.hh"
#include "NodeList/nthNodalMoment.hh"
#include "Kernel/TableKernel.hh"
#include "Mesh/Mesh.hh"

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
namespace Spheral {

typedef NodeListRegistrar<Dim<1> > NodeListRegistrar1d;
typedef NodeListRegistrar<Dim<2> > NodeListRegistrar2d;
typedef NodeListRegistrar<Dim<3> > NodeListRegistrar3d;

typedef NodeList<Dim<1> > NodeList1d;
typedef NodeList<Dim<2> > NodeList2d;
typedef NodeList<Dim<3> > NodeList3d;

typedef FluidNodeList<Dim<1> > FluidNodeList1d;
typedef FluidNodeList<Dim<2> > FluidNodeList2d;
typedef FluidNodeList<Dim<3> > FluidNodeList3d;

typedef SolidNodeList<Dim<1> > SolidNodeList1d;
typedef SolidNodeList<Dim<2> > SolidNodeList2d;
typedef SolidNodeList<Dim<3> > SolidNodeList3d;

typedef SmoothingScaleBase<Dim<1> > SmoothingScaleBase1d;
typedef SmoothingScaleBase<Dim<2> > SmoothingScaleBase2d;
typedef SmoothingScaleBase<Dim<3> > SmoothingScaleBase3d;

typedef FixedSmoothingScale<Dim<1> > FixedSmoothingScale1d;
typedef FixedSmoothingScale<Dim<2> > FixedSmoothingScale2d;
typedef FixedSmoothingScale<Dim<3> > FixedSmoothingScale3d;

typedef SPHSmoothingScale<Dim<1> > SPHSmoothingScale1d;
typedef SPHSmoothingScale<Dim<2> > SPHSmoothingScale2d;
typedef SPHSmoothingScale<Dim<3> > SPHSmoothingScale3d;

typedef ASPHSmoothingScale<Dim<1> > ASPHSmoothingScale1d;
typedef ASPHSmoothingScale<Dim<2> > ASPHSmoothingScale2d;
typedef ASPHSmoothingScale<Dim<3> > ASPHSmoothingScale3d;

// //------------------------------------------------------------------------------
// // Provide a non-iterator based interface to generateVoidNodes.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// generateVoidNodes(const std::vector<NodeList<Dimension>*>& nodeLists,
//                   const Mesh<Dimension>& mesh,
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
FieldList<Dimension, typename MomentTraits<Dimension, moment>::Moment>
nthNodalMoment(const std::vector<NodeList<Dimension>*>& nodeLists,
               const TableKernel<Dimension>& W,
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
                           const TableKernel<Dimension>& W,
                           const bool useGradientAsKernel,
                           FieldList<Dimension, typename Dimension::Scalar>& zerothMoment,
                           FieldList<Dimension, typename Dimension::Vector>& firstMoment) {
  zerothAndFirstNodalMoments<Dimension, typename std::vector<NodeList<Dimension>*>::const_iterator>
    (nodeLists.begin(), nodeLists.end(), W, useGradientAsKernel, zerothMoment, firstMoment);
}

}

typedef std::vector<Spheral::NodeList<Spheral::Dim<1> >*> vector_of_NodeList1d;
typedef std::vector<Spheral::NodeList<Spheral::Dim<2> >*> vector_of_NodeList2d;
typedef std::vector<Spheral::NodeList<Spheral::Dim<3> >*> vector_of_NodeList3d;

typedef std::vector<const Spheral::NodeList<Spheral::Dim<1> >*> vector_of_const_NodeList1d;
typedef std::vector<const Spheral::NodeList<Spheral::Dim<2> >*> vector_of_const_NodeList2d;
typedef std::vector<const Spheral::NodeList<Spheral::Dim<3> >*> vector_of_const_NodeList3d;

typedef std::vector<Spheral::FluidNodeList<Spheral::Dim<1> >*> vector_of_FluidNodeList1d;
typedef std::vector<Spheral::FluidNodeList<Spheral::Dim<2> >*> vector_of_FluidNodeList2d;
typedef std::vector<Spheral::FluidNodeList<Spheral::Dim<3> >*> vector_of_FluidNodeList3d;

typedef std::vector<const Spheral::FluidNodeList<Spheral::Dim<1> >*> vector_of_const_FluidNodeList1d;
typedef std::vector<const Spheral::FluidNodeList<Spheral::Dim<2> >*> vector_of_const_FluidNodeList2d;
typedef std::vector<const Spheral::FluidNodeList<Spheral::Dim<3> >*> vector_of_const_FluidNodeList3d;

typedef std::pair<const Spheral::NodeList<Spheral::Dim<1> >*, std::string> pair_NodeList1d_string;
typedef std::pair<const Spheral::NodeList<Spheral::Dim<2> >*, std::string> pair_NodeList2d_string;
typedef std::pair<const Spheral::NodeList<Spheral::Dim<3> >*, std::string> pair_NodeList3d_string;

typedef std::vector<pair_NodeList1d_string> vector_of_pair_NodeList1d_string;
typedef std::vector<pair_NodeList2d_string> vector_of_pair_NodeList2d_string;
typedef std::vector<pair_NodeList3d_string> vector_of_pair_NodeList3d_string;

typedef std::vector<Spheral::NodeList<Spheral::Dim<1> >*>::iterator vector_of_NodeList1d_iterator;
typedef std::vector<Spheral::NodeList<Spheral::Dim<2> >*>::iterator vector_of_NodeList2d_iterator;
typedef std::vector<Spheral::NodeList<Spheral::Dim<3> >*>::iterator vector_of_NodeList3d_iterator;

typedef std::vector<Spheral::FluidNodeList<Spheral::Dim<1> >*>::iterator vector_of_FluidNodeList1d_iterator;
typedef std::vector<Spheral::FluidNodeList<Spheral::Dim<2> >*>::iterator vector_of_FluidNodeList2d_iterator;
typedef std::vector<Spheral::FluidNodeList<Spheral::Dim<3> >*>::iterator vector_of_FluidNodeList3d_iterator;

typedef std::vector<Spheral::SolidNodeList<Dim<1> >*> vector_of_SolidNodeList1d;
typedef std::vector<Spheral::SolidNodeList<Dim<2> >*> vector_of_SolidNodeList2d;
typedef std::vector<Spheral::SolidNodeList<Dim<3> >*> vector_of_SolidNodeList3d;

typedef std::vector<const Spheral::SolidNodeList<Dim<1> >*> vector_of_const_SolidNodeList1d;
typedef std::vector<const Spheral::SolidNodeList<Dim<2> >*> vector_of_const_SolidNodeList2d;
typedef std::vector<const Spheral::SolidNodeList<Dim<3> >*> vector_of_const_SolidNodeList3d;

typedef std::vector<Spheral::SolidNodeList<Dim<1> >*>::iterator vector_of_SolidNodeList1d_iterator;
typedef std::vector<Spheral::SolidNodeList<Dim<2> >*>::iterator vector_of_SolidNodeList2d_iterator;
typedef std::vector<Spheral::SolidNodeList<Dim<3> >*>::iterator vector_of_SolidNodeList3d_iterator;

//------------------------------------------------------------------------------
// Extract the NodeListRegistrar instance.
//------------------------------------------------------------------------------
namespace Spheral {

template<typename Dimension>
inline
static NodeListRegistrar<Dimension>*
nodeListRegistrarInstance() {
  return &(NodeListRegistrar<Dimension>::instance());
}

}

#endif
