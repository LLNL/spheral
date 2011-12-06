//---------------------------------Spheral++----------------------------------//
// Mesh -- the generic mesh class we construct around meshless nodes.
// 1D:  Mesh is a collection of line segments.
// 2D:  Arbitrary polygons.
// 3D:  Arbitrary polyhedra.
//
// Created by JMO, Thu Oct 14 11:08:47 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_Mesh__
#define __Spheral_Mesh__

#include <stdint.h>
#include <vector>
#include <map>
#include <string>
#include "boost/shared_ptr.hpp"
#include "boost/tuple/tuple.hpp"

#include "Geometry/Dimension.hh"
#include "VoroPP.hh"
#include "MeshWall.hh"

namespace Spheral {

namespace NodeSpace {
  template<typename Dimension> class NodeList;
}

namespace MeshSpace {

template<typename Dimension>
class Mesh {
  //--------------------------- Public Interface ---------------------------//
public:
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::ConvexHull ConvexHull;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef uint64_t KeyElement;
  typedef boost::tuple<KeyElement, KeyElement, KeyElement> Key;
  static const unsigned UNSETID;
  static const unsigned minFacesPerZone, minEdgesPerZone, minNodesPerZone;
  static const unsigned minEdgesPerFace, minNodesPerFace;

  //---------------------------------------------------------------------------
  // Declare the elements of the mesh.
  //---------------------------------------------------------------------------
  class Node;
  class Edge;
  class Face;
  class Zone;

  //---------------------------------------------------------------------------
  // Types.
  //---------------------------------------------------------------------------
  typedef std::vector<Node> NodeContainer;
  typedef std::vector<Edge> EdgeContainer;
  typedef std::vector<Face> FaceContainer;
  typedef std::vector<Zone> ZoneContainer;

  typedef typename NodeContainer::const_iterator NodeIterator;
  typedef typename EdgeContainer::const_iterator EdgeIterator;
  typedef typename FaceContainer::const_iterator FaceIterator;
  typedef typename ZoneContainer::const_iterator ZoneIterator;

  typedef boost::shared_ptr<MeshWall<Dimension> > MeshWallPtr;

  //---------------------------------------------------------------------------
  // Static methods.
  //---------------------------------------------------------------------------
  static const int nDim() { return Dimension::nDim; }
  
  //---------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //---------------------------------------------------------------------------
  Mesh();
  Mesh(const std::vector<Vector>& generators,
       const Vector& xmin,
       const Vector& xmax);
  Mesh(const std::vector<Vector>& generators,
       const FacetedVolume& boundary);
  Mesh(const std::vector<Vector>& nodePositions,
       const std::vector<std::vector<unsigned> >& edgeNodes,
       const std::vector<std::vector<unsigned> >& faceEdges,
       const std::vector<std::vector<unsigned> >& zoneFaces);
  Mesh& operator=(const Mesh& rhs);
  ~Mesh();

  // Clear out any exising data in the mesh.
  void clear();

  // Add walls to this mesh (for generation).
  void addWall(MeshWallPtr wallPtr);

  // Reconstruct from a set of generators.
  void reconstruct(const std::vector<Vector>& generators,
                   const Vector& xmin,
                   const Vector& xmax);
  void reconstruct(const std::vector<Vector>& generators,
                   const FacetedVolume& boundary);
  template<typename BoundaryIterator>
  void reconstruct(const std::vector<Vector>& generators,
                   const Vector& xmin,
                   const Vector& xmax,
                   const BoundaryIterator boundaryBegin,
                   const BoundaryIterator boundaryEnd);
  template<typename BoundaryIterator>
  void reconstruct(const std::vector<Vector>& generators,
                   const FacetedVolume& boundary,
                   const BoundaryIterator boundaryBegin,
                   const BoundaryIterator boundaryEnd);

  // Remove zones from the mesh according to a mask:
  //   mask[i] = 0  ---> remove zone i
  //   mask[i] = 1  ---> keep zone i
  void removeZonesByMask(const std::vector<unsigned>& mask);

  // Remove edges below a threshold fraction size.
  void cleanEdges(const double edgeTol);

  // Sizes.
  unsigned numNodes() const;
  unsigned numEdges() const;
  unsigned numFaces() const;
  unsigned numZones() const;

  // Access the elements.
  const Node& node(const unsigned i) const;
  NodeIterator nodeBegin() const;
  NodeIterator nodeEnd() const;

  const Edge& edge(const unsigned i) const;
  EdgeIterator edgeBegin() const;
  EdgeIterator edgeEnd() const;

  const Face& face(const unsigned i) const;
  FaceIterator faceBegin() const;
  FaceIterator faceEnd() const;

  const Zone& zone(const unsigned i) const;
  ZoneIterator zoneBegin() const;
  ZoneIterator zoneEnd() const;

  // We also provide the ability to extract the zone corresponding to the given node
  // in a NodeList.
  const Zone& zone(const NodeSpace::NodeList<Dimension>& nodeList, const unsigned i) const;
  const Zone& zone(const unsigned nodeListi, const unsigned i) const;

  // Extract the zone offset for the given NodeList.
  unsigned offset(const NodeSpace::NodeList<Dimension>& nodeList) const;
  unsigned offset(const unsigned nodeListi) const;

  // Compute the parallel connectivity.
  void generateDomainInfo();

  // Check that the internal parallel info is consistent.
  std::string validDomainInfo(const Vector& xmin,
                              const Vector& xmax,
                              const bool checkUniqueSendProc) const;

  // Access to the parallel info.
  const std::vector<unsigned>& neighborDomains() const;
  const std::vector<std::vector<unsigned> >& sharedNodes() const;

  // Compute the minimum scale (distance between nodes).
  double minimumScale() const;

  // Store the given offsets for a set of NodeLists.
  template<typename NodeListIterator>
  void storeNodeListOffsets(const NodeListIterator nodeListBegin, 
                            const NodeListIterator nodeListEnd,
                            const std::vector<unsigned>& offsets);
  void storeNodeListOffsets(const std::vector<NodeSpace::NodeList<Dimension>*>& nodeListPtrs,
                            const std::vector<unsigned>& offsets);

  //--------------------------- Private Interface ---------------------------//
private:
  // The mesh data.
  std::vector<Vector> mNodePositions;
  NodeContainer mNodes;
  EdgeContainer mEdges;
  FaceContainer mFaces;
  ZoneContainer mZones;

  // The optional parallel info.
  std::vector<unsigned> mNeighborDomains;
  std::vector<std::vector<unsigned> > mSharedNodes;

  // The offsets into the zones for each NodeList.
  std::map<std::string, unsigned> mNodeListNameOffsets;
  std::vector<unsigned> mNodeListIndexOffsets;

  // The set of MeshWalls.
  std::vector<MeshWallPtr> mWallPtrs;

  // Internal method to recompute IDs based on a mask indicating which 
  // elements are to be kept/removed.
  std::vector<unsigned> recomputeIDs(const std::vector<unsigned>& mask) const;

  // Reassign IDs in a vector of IDs based on an map of old -> new IDs.
  void reassignIDs(std::vector<unsigned>& ids,
                   const std::vector<unsigned>& old2new) const;

  // Internal method to handle reconstructing the mesh after boundary
  // conditions and such have all been provided.
  void reconstructInternal(const std::vector<Vector>& generators,
                           const Vector& xmin, 
                           const Vector& xmax);
};

}
}

#include "Node.hh"
#include "Edge.hh"
#include "Face.hh"
#include "Zone.hh"
#include "MeshInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace MeshSpace {
    template<typename Dimension> class Mesh;
  }
}

#endif
