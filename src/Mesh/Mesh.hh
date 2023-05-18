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
#include <tuple>
#include <memory>

#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension> class NodeList;

template<typename Dimension>
class Mesh {
  //--------------------------- Public Interface ---------------------------//
public:
  typedef Dimension DimensionType;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ConvexHull ConvexHull;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef uint64_t KeyElement;
  typedef std::tuple<KeyElement, KeyElement, KeyElement> Key;
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

  //---------------------------------------------------------------------------
  // Static methods.
  //---------------------------------------------------------------------------
  static int nDim() { return Dimension::nDim; }
  
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
       const std::vector<std::vector<int> >& zoneFaces);
  Mesh(const Mesh& rhs);
  Mesh& operator=(const Mesh& rhs);
  virtual ~Mesh();

  // Clear out any exising data in the mesh.
  void clear();

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
  const Face& face(const int i) const;
  FaceIterator faceBegin() const;
  FaceIterator faceEnd() const;

  const Zone& zone(const unsigned i) const;
  const Zone& zone(const int i) const;
  ZoneIterator zoneBegin() const;
  ZoneIterator zoneEnd() const;

  // We also provide the ability to extract the zone corresponding to the given node
  // in a NodeList.
  const Zone& zone(const NodeList<Dimension>& nodeList, const unsigned i) const;
  const Zone& zone(const unsigned nodeListi, const unsigned i) const;

  // Extract the zone offset for the given NodeList.
  unsigned offset(const NodeList<Dimension>& nodeList) const;
  unsigned offset(const unsigned nodeListi) const;

  // Look up the (nodeListID, nodeID) corresponding to the given zoneID.
  void lookupNodeListID(const unsigned zoneID, unsigned& nodeListi, unsigned& i) const;

  // Compute the communicated mesh structures.
  void generateDomainInfo();

  // Generate a parallel rind of cells around each domain representing a one zone
  // thick set of zones shared with the neighboring processors.
  // Note we do not recompute the shared elements (nodes & faces) as part of this
  // procedure, so following this operation those shared elements are no longer
  // on the surface of the local mesh!
  void generateParallelRind();

  // This version also exchanges the generators for the rind cells.
  void generateParallelRind(std::vector<Vector>& /*generators*/,
                            std::vector<SymTensor>& /*Hs*/);

  // Compute unique global IDs for each node.
  std::vector<unsigned> globalMeshNodeIDs() const;

  // Compute unique global IDs for each face.
  std::vector<unsigned> globalMeshFaceIDs(const std::vector<unsigned>& globalNodeIDs) const;

  // Access to the parallel info.
  const std::vector<unsigned>& neighborDomains() const;
  const std::vector<std::vector<unsigned> >& sharedNodes() const;
  const std::vector<std::vector<unsigned> >& sharedFaces() const;

  // Compute the minimum scale (distance between nodes).
  double minimumScale() const;

  // Store the given offsets for a set of NodeLists.
  template<typename NodeListIterator>
  void storeNodeListOffsets(const NodeListIterator nodeListBegin, 
                            const NodeListIterator nodeListEnd,
                            const std::vector<unsigned>& offsets);
  void storeNodeListOffsets(const std::vector<NodeList<Dimension>*>& nodeListPtrs,
                            const std::vector<unsigned>& offsets);

  // Compute the bounding surface of the mesh.
  FacetedVolume boundingSurface() const;
  void boundingBox(Vector& xmin, Vector& xmax) const;

  // Encapsulate the ones complement for signed (oriented) IDs.
  static unsigned int positiveID(const int id);

  // Perform basic mesh validity checks.
  virtual std::string valid() const;

  // Check that the internal parallel info is consistent.
  virtual std::string validDomainInfo(const Vector& xmin,
                                      const Vector& xmax,
                                      const bool /*checkUniqueSendProc*/) const;

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
  std::vector<std::vector<unsigned> > mSharedNodes, mSharedFaces;

  // The offsets into the zones for each NodeList.
  std::map<std::string, unsigned> mNodeListNameOffsets;
  std::vector<unsigned> mNodeListIndexOffsets;

  // Internal method to recompute IDs based on a mask indicating which 
  // elements are to be kept/removed.
  std::vector<unsigned> recomputeIDs(const std::vector<unsigned>& mask) const;

  // Reassign IDs in a vector of IDs based on an map of old -> new IDs.
  void reassignIDs(std::vector<unsigned>& ids,
                   const std::vector<unsigned>& old2new) const;
  void reassignIDs(std::vector<int>& ids,
                   const std::vector<unsigned>& old2new) const;

  // Delete any elements in the list which are set to UNSETID.
  void removeUNSETIDs(std::vector<unsigned>& ids) const;
  void removeUNSETIDs(std::vector<int>& ids) const;

  // Internal methods to handle reconstructing the mesh after boundary
  // conditions and such have all been provided.
  void reconstructInternal(const std::vector<Vector>& generators,
                           const Vector& xmin, 
                           const Vector& xmax);
  void reconstructInternal(const std::vector<Vector>& generators,
                           const FacetedVolume& boundary);

  // Internal method used to add on new mesh elements based on sets of 
  // nodes arranged into Face arrays.
  void createNewMeshElements(const std::vector<std::vector<std::vector<unsigned> > >& newCells);
};

// Declare 1D specializations.
template<> inline void Mesh<Dim<1> >::cleanEdges(const double /*edgeTol*/) {}

}

#include "Node.hh"
#include "Edge.hh"
#include "Face.hh"
#include "Zone.hh"
#include "MeshInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class Mesh;
}

#endif
