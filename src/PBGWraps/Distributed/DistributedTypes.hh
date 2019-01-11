#ifndef __PBGWRAPS_DISTRIBUTEDTYPES__
#define __PBGWRAPS_DISTRIBUTEDTYPES__

#include "referenceAsPointer.hh"
#include "Distributed/DistributedBoundary.hh"
#include "Distributed/NestedGridDistributedBoundary.hh"
#include "Distributed/BoundingVolumeDistributedBoundary.hh"
#include "Distributed/TreeDistributedBoundary.hh"
#include "Distributed/DomainNode.hh"
#include "Distributed/RedistributeNodes.hh"
#include "Distributed/DistributeByXPosition.hh"
// #include "Distributed/NestedGridRedistributeNodes.hh"
#include "Distributed/SpaceFillingCurveRedistributeNodes.hh"
#include "Distributed/MortonOrderRedistributeNodes.hh"
#include "Distributed/PeanoHilbertOrderRedistributeNodes.hh"
#include "Distributed/SortAndDivideRedistributeNodes1d.hh"
#include "Distributed/SortAndDivideRedistributeNodes2d.hh"
#include "Distributed/SortAndDivideRedistributeNodes3d.hh"
#include "Distributed/VoronoiRedistributeNodes.hh"

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
namespace Spheral {

typedef DistributedBoundary<Dim<1> > DistributedBoundary1d;
typedef DistributedBoundary<Dim<2> > DistributedBoundary2d;
typedef DistributedBoundary<Dim<3> > DistributedBoundary3d;

typedef DistributedBoundary<Dim<1> >::DomainBoundaryNodes DomainBoundaryNodes1d;
typedef DistributedBoundary<Dim<2> >::DomainBoundaryNodes DomainBoundaryNodes2d;
typedef DistributedBoundary<Dim<3> >::DomainBoundaryNodes DomainBoundaryNodes3d;

typedef NestedGridDistributedBoundary<Dim<1> > NestedGridDistributedBoundary1d;
typedef NestedGridDistributedBoundary<Dim<2> > NestedGridDistributedBoundary2d;
typedef NestedGridDistributedBoundary<Dim<3> > NestedGridDistributedBoundary3d;

typedef BoundingVolumeDistributedBoundary<Dim<1> > BoundingVolumeDistributedBoundary1d;
typedef BoundingVolumeDistributedBoundary<Dim<2> > BoundingVolumeDistributedBoundary2d;
typedef BoundingVolumeDistributedBoundary<Dim<3> > BoundingVolumeDistributedBoundary3d;

typedef TreeDistributedBoundary<Dim<1> > TreeDistributedBoundary1d;
typedef TreeDistributedBoundary<Dim<2> > TreeDistributedBoundary2d;
typedef TreeDistributedBoundary<Dim<3> > TreeDistributedBoundary3d;

typedef DomainNode<Dim<1> > DomainNode1d;
typedef DomainNode<Dim<2> > DomainNode2d;
typedef DomainNode<Dim<3> > DomainNode3d;

typedef RedistributeNodes<Dim<1> > RedistributeNodes1d;
typedef RedistributeNodes<Dim<2> > RedistributeNodes2d;
typedef RedistributeNodes<Dim<3> > RedistributeNodes3d;

typedef DistributeByXPosition<Dim<1> > DistributeByXPosition1d;
typedef DistributeByXPosition<Dim<2> > DistributeByXPosition2d;

// typedef NestedGridRedistributeNodes<Dim<1> > NestedGridRedistributeNodes1d;
// typedef NestedGridRedistributeNodes<Dim<2> > NestedGridRedistributeNodes2d;
// typedef NestedGridRedistributeNodes<Dim<3> > NestedGridRedistributeNodes3d;

typedef SpaceFillingCurveRedistributeNodes<Dim<1> > SpaceFillingCurveRedistributeNodes1d;
typedef SpaceFillingCurveRedistributeNodes<Dim<2> > SpaceFillingCurveRedistributeNodes2d;
typedef SpaceFillingCurveRedistributeNodes<Dim<3> > SpaceFillingCurveRedistributeNodes3d;

typedef MortonOrderRedistributeNodes<Dim<1> > MortonOrderRedistributeNodes1d;
typedef MortonOrderRedistributeNodes<Dim<2> > MortonOrderRedistributeNodes2d;
typedef MortonOrderRedistributeNodes<Dim<3> > MortonOrderRedistributeNodes3d;

typedef PeanoHilbertOrderRedistributeNodes<Dim<1> > PeanoHilbertOrderRedistributeNodes1d;
typedef PeanoHilbertOrderRedistributeNodes<Dim<2> > PeanoHilbertOrderRedistributeNodes2d;
typedef PeanoHilbertOrderRedistributeNodes<Dim<3> > PeanoHilbertOrderRedistributeNodes3d;

typedef VoronoiRedistributeNodes<Dim<1> > VoronoiRedistributeNodes1d;
typedef VoronoiRedistributeNodes<Dim<2> > VoronoiRedistributeNodes2d;
typedef VoronoiRedistributeNodes<Dim<3> > VoronoiRedistributeNodes3d;

typedef std::pair<uint64_t, DomainNode1d> pair_ULL_DomainNode1d;
typedef std::pair<uint64_t, DomainNode2d> pair_ULL_DomainNode2d;
typedef std::pair<uint64_t, DomainNode3d> pair_ULL_DomainNode3d;

typedef std::vector<DomainNode1d> vector_of_DomainNode1d;
typedef std::vector<DomainNode2d> vector_of_DomainNode2d;
typedef std::vector<DomainNode3d> vector_of_DomainNode3d;

typedef std::vector<pair_ULL_DomainNode1d> vector_of_pair_ULL_DomainNode1d;
typedef std::vector<pair_ULL_DomainNode2d> vector_of_pair_ULL_DomainNode2d;
typedef std::vector<pair_ULL_DomainNode3d> vector_of_pair_ULL_DomainNode3d;

}

// namespace Spheral {

// //------------------------------------------------------------------------------
// // Get the NestedGridNeighbor from a NodeList.
// //------------------------------------------------------------------------------

// template<typename Dimension>
// inline
// NestedGridNeighbor<Dimension>*
// getNestedGridNeighbor(const DistributedBoundary<Dimension>& self,
//                       const NodeList<Dimension>& nodeList) {
//   return &(self.getNestedGridNeighbor(&nodeList));
// }

// //------------------------------------------------------------------------------
// // Get the NestedGridDistributedBoundary instance.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// static
// NestedGridDistributedBoundary<Dimension>*
// getNestedGridDistributedBoundaryInstance() {
//   return &(NestedGridDistributedBoundary<Dimension>::instance());
// }

// //------------------------------------------------------------------------------
// // Get the TreeDistributedBoundary instance.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// static
// NestedGridDistributedBoundary<Dimension>*
// getTreeDistributedBoundaryInstance() {
//   return &(TreeDistributedBoundary<Dimension>::instance());
// }

// //------------------------------------------------------------------------------
// // Get the BoundingVolumeDistributedBoundary instance.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// static
// BoundingVolumeDistributedBoundary<Dimension>*
// getBoundingVolumeDistributedBoundaryInstance() {
//   return &(BoundingVolumeDistributedBoundary<Dimension>::instance());
// }

// }

#endif
