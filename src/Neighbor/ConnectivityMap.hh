//---------------------------------Spheral++----------------------------------//
// ConnectivityMap
//
// Stores the full set of significant neighbors for a set of NodeLists.
//
// Created by J. Michael Owen, Sun Oct 30 15:36:33 PST 2005
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_ConnectivityMap_hh_
#define _Spheral_NeighborSpace_ConnectivityMap_hh_

#ifndef __GCCXML__
#include <vector>
#include <map>
#include "Utilities/KeyTraits.hh"
#include "Field/FieldList.hh"
#else
#include "fakestl.hh"
#endif

namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
    template<typename Dimension> class FluidNodeList;
  }
  namespace BoundarySpace {
    template<typename Dimension> class Boundary;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }
}

namespace Spheral {
namespace NeighborSpace {

template<typename Dimension>
class ConnectivityMap {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef std::vector<int>::const_iterator const_iterator;

  // Constructors, destructor.
  ConnectivityMap();
  ~ConnectivityMap();

  template<typename NodeListIterator>
  ConnectivityMap(const NodeListIterator& begin,
                  const NodeListIterator& end,
                  const bool buildGhostConnectivity);

  // Rebuild for a given set of NodeLists.
  template<typename NodeListIterator>
  void rebuild(const NodeListIterator& begin, 
               const NodeListIterator& end, 
               const bool computeGhostConnectivity);

  // Patch the connectivity information:
  // flags   -- (0,1): 0 => node deleted, 1 => node preserved
  // old2new -- maps old -> new node indices.
  void patchConnectivity(const FieldSpace::FieldList<Dimension, int>& flags,
                         const FieldSpace::FieldList<Dimension, int>& old2new);

  // Are we computing neighbors for ghosts?
  bool buildGhostConnectivity() const;

  // Get the set of NodeLists.
  const std::vector<const NodeSpace::NodeList<Dimension>*>& nodeLists() const;

  // Get the set of neighbors for the given (internal!) node in the given NodeList.
  const std::vector< std::vector<int> >&
  connectivityForNode(const NodeSpace::NodeList<Dimension>* nodeListPtr,
                      const int nodeID) const;

  // Same as above, just referencing the NodeList by an integer index.
  const std::vector< std::vector<int> >&
  connectivityForNode(const int nodeListID,
                      const int nodeID) const;

  // Compute the common neighbors for a pair of nodes.  Note this method 
  // returns by value since this information is not stored by ConnectivityMap.
  std::vector< std::vector<int> >
  connectivityIntersectionForNodes(const int nodeListi, const int i,
                                   const int nodeListj, const int j) const;

  // Compute the union of neighbors for a pair of nodes.  Note this method 
  // returns by value since this information is not stored by ConnectivityMap.
  std::vector< std::vector<int> >
  connectivityUnionForNodes(const int nodeListi, const int i,
                            const int nodeListj, const int j) const;

  // Compute the number of neighbors for the given node.
  int numNeighborsForNode(const NodeSpace::NodeList<Dimension>* nodeListPtr,
                          const int nodeID) const;

  int numNeighborsForNode(const int nodeListID,
                          const int nodeID) const;

  // Return the connectivity in terms of global node IDs.
  std::map<int, std::vector<int> > globalConnectivity(std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries) const;

  // Function to determine if given node information (i and j), if the 
  // pair should already have been calculated by iterating over each
  // others neighbors.
#pragma omp declare target
  bool calculatePairInteraction(const int nodeListi, const int i, 
                                const int nodeListj, const int j,
                                const int firstGhostNodej) const;

#pragma omp end declare target
  // Provide iterator interface for walking the nodes in a NodeList
  // in order to maintain domain decomposition independence when 
  // desired.
  const_iterator begin(const int nodeList) const;
  const_iterator end(const int nodeList) const;

  // Get the ith NodeList or FluidNodeList.
  const NodeSpace::NodeList<Dimension>& nodeList(const int index) const;

  // Return which NodeList index in order the given one would be in our connectivity.
  unsigned nodeListIndex(const NodeSpace::NodeList<Dimension>* nodeListPtr) const;

  // Check that the internal data structure is valid.
  bool valid() const;

private:
#ifndef __GCCXML__
  //--------------------------- Private Interface ---------------------------//
  // The set of NodeLists.
  std::vector<const NodeSpace::NodeList<Dimension>*> mNodeLists;

  // Are we building ghost connectivity?
  bool mBuildGhostConnectivity;

  // The full connectivity map.  This might be quite large!
  // [offset[NodeList] + nodeID] [NodeListID] [neighborIndex]
  typedef std::vector<std::vector<std::vector<int> > > ConnectivityStorageType;
  std::vector<int> mOffsets;
  ConnectivityStorageType mConnectivity;

  // The set of node indices per Nodelist in order for traversal.
  std::vector< std::vector<int> > mNodeTraversalIndices;

  // The set of keys we may compute for each node.
  typedef typename KeyTraits::Key Key;
  FieldSpace::FieldList<Dimension, Key> mKeys;

  // Internal method to fill in the connectivity, once the set of NodeLists 
  // is determined.
  void computeConnectivity();

  // No default constructor, copying, or assignment.
  ConnectivityMap(const ConnectivityMap&);
  ConnectivityMap& operator=(const ConnectivityMap&);
#endif
};

}
}

#ifndef __GCCXML__
#include "ConnectivityMapInline.hh"
#endif

#endif

