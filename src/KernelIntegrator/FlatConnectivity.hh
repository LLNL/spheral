//---------------------------------Spheral++----------------------------------//
// FlatConnectivity
//
// Creates global indices for each point and flattened connectivity indices
//----------------------------------------------------------------------------//
#ifndef __Spheral_FlatConnectivity_hh__
#define __Spheral_FlatConnectivity_hh__

#include <vector>
#include <unordered_map>

#include "DataBase/DataBase.hh"

namespace Spheral {

// These hashes are specialized for low numbers of integers
// The array version assumes the last index is either 0 or 1
template<typename DataType> struct NormalHash {
  int operator()(const DataType& v) const;
};

template<typename Dimension>
class FlatConnectivity {
public:
  // Return this if connectivity does not exist
  static constexpr int NoConnectivity = -1;
  
  // Typedefs
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef typename std::unordered_map<int, int> IndexMap;
  typedef typename std::array<int, Dimension::nDim> ArrayDim;
  typedef typename std::unordered_map<ArrayDim, int, NormalHash<ArrayDim>> NormalMap;

  // Constructor
  FlatConnectivity();

  // Check whether things are initialized
  bool indexingInitialized() const;
  bool overlapIndexingInitialized() const;
  bool globalIndexingInitialized() const;
  bool surfaceIndexingInitialized() const;
  bool boundaryInformationInitialized() const;
  
  // Lower and upper global indices
  int firstGlobalIndex() const;
  int lastGlobalIndex() const;
  
  // Local and global number of elements
  int numNodes() const;
  int numInternalNodes() const;
  int numGlobalNodes() const;
  int numBoundaryNodes() const;
  
  // Get local index from NodeList and Node indices
  int nodeToLocal(const int nodeListi, const int nodei) const;
  
  // Get NodeList and Node indices from local index
  std::pair<int, int> localToNode(const int locali) const;

  // Get global index from local index
  int localToGlobal(const int locali) const;
  
  // Number of neighbors, including self, for this point
  int numNeighbors(const int locali) const;
  int numOverlapNeighbors(const int locali) const;
  int numConstNeighbors(const int locali) const;
  int numConstOverlapNeighbors(const int locali) const;
  int numNonConstNeighbors(const int locali) const;
  int numNonConstOverlapNeighbors(const int locali) const;
  
  // For the point i, for its neighbor j, get the flattened index for point j
  int localToFlat(const int locali, const int localj) const; // return flatj
  int localToFlatOverlap(const int locali, const int localj) const; // return flatj

  // For the point i, for the flattened index for j, get the local index for j
  int flatToLocal(const int locali, const int flatj) const;  // return localj
  int flatOverlapToLocal(const int locali, const int flatj) const; // return localj

  // Is the point i a constant boundary node?
  bool isConstantBoundaryNode(const int locali) const;
  
  // Get the local indices for the neighbors of point i
  void neighborIndices(const int locali,
                       std::vector<int>& localNeighbors) const;
  void overlapNeighborIndices(const int locali,
                              std::vector<int>& localNeighbors) const;
  void constNeighborIndices(const int locali,
                            std::vector<int>& localNeighbors) const;
  void overlapConstNeighborIndices(const int locali,
                                   std::vector<int>& localNeighbors) const;
  void nonConstNeighborIndices(const int locali,
                               std::vector<int>& localNeighbors) const;
  void overlapNonConstNeighborIndices(const int locali,
                                      std::vector<int>& localNeighbors) const;
  
  // Get the global indices for the neighbors of point i, not including any constant boundary nodes
  void globalNeighborIndices(const int locali,
                             std::vector<int>& globalNeighborIndices) const;
  void globalOverlapNeighborIndices(const int locali,
                                    std::vector<int>& globalNeighborIndices) const;
  
  // Get the number of unique void surfaces included in this point's neighbors
  int numSurfaces(const int locali) const;
  
  // Get the unique void surface index from the normal direction
  int surfaceIndex(const int locali,
                   const Vector& normal) const; // return flats
  
  // Same as above, but using the key for better efficiency
  int surfaceIndex(const int locali,
                   const ArrayDim& values) const; // return flats
  
  // For the point i, get the normal directions for each void surface
  const Vector& normal(const int locali, const int flats) const;

  // For the voronoi cell i, get the number of void surfaces
  int numSurfacesForCell(const int locali) const;
  
  // For the voronoi cell i, get the surface index for this void surface
  int surfaceIndexForCell(const int locali, const int flats) const;

  // Convert a normal to an array of ints and back, for surface indexing
  void normalToArray(const Vector& normal,
                     ArrayDim& values) const;
  void arrayToNormal(const ArrayDim& values,
                     Vector& normal) const;
  
  // Compute the local point indices:
  // call this first
  void computeIndices(const DataBase<Dimension>& dataBase);

  // Compute the local overlap connectivity:
  // call after computeIndices
  void computeOverlapIndices(const DataBase<Dimension>& dataBase);

  // Compute the global point indices, optionally applying the boundaries to the global indices:
  // call after computeIndices
  void computeGlobalIndices(const DataBase<Dimension>& dataBase,
                            const std::vector<Boundary<Dimension>*>& boundaries); 
  
  // Compute the surface indices:
  // call after computeIndices
  void computeSurfaceIndices(const DataBase<Dimension>& dataBase,
                             const State<Dimension>& state);
  
  // Compute the constant boundary nodes: 
  // call after computeIndices and computeOverlapIndices
  void computeBoundaryInformation(const DataBase<Dimension>& dataBase,
                                  const std::vector<Boundary<Dimension>*>& boundaries);

  // // Check whether two functions overlap
  // bool checkOverlap(const Vector& x1,
  //                   const SymTensor& H1,
  //                   const Vector& x2,
  //                   const SymTensor& H2,
  //                   const Scalar extent) const;
  
private:
  
  // Keep track of what has been initialized
  bool mIndexingInitialized;
  bool mGhostIndexingInitialized;
  bool mOverlapIndexingInitialized;
  bool mGlobalIndexingInitialized;
  bool mSurfaceIndexingInitialized;
  bool mBoundaryInformationInitialized;

  // Size information
  int mNumLocalNodes;
  int mNumInternalLocalNodes;
  int mNumConnectivityNodes;
  int mNumGlobalNodes;
  int mNumBoundaryNodes;
  int mFirstGlobalIndex;
  int mLastGlobalIndex;

  /* To make the indices clearer, they are commented in terms of these variables:
     nodeListi, nodei: The standard node indices
     locali: The flattened node index for this processor
     globali: The flattened node index globally
     flatj: The flattened index for the jth neighbor of locali
     flats: The flattened surface index s for point i
     normal: The normal direction for the surface
     normalasarray: The normal direction, converted to an array for indexing */
  
  // Flat indexing: a single index for each local point
  std::vector<std::vector<int>> mNodeToLocalIndex; // (nodeListi, nodei) -> locali
  std::vector<std::pair<int, int>> mLocalToNodeIndex; // locali -> (nodeListi, nodei)

  // Global indexing
  std::vector<int> mLocalToGlobalIndex; // locali -> globali
  
  // Connectivity: points that see each other's center
  std::vector<int> mNumNeighbors; // locali -> numNeighborsi
  std::vector<std::vector<int>> mFlatToLocalIndex; // (locali, flatj) -> localj
  std::vector<IndexMap> mLocalToFlatIndex; // (locali, localj) -> flatj

  // Overlap connectivity: points whose kernels overlap
  std::vector<int> mNumOverlapNeighbors; // localli -> numOverlapNeighborsi
  std::vector<std::vector<int>> mFlatOverlapToLocalIndex; // (locali, flatj) -> localj
  std::vector<IndexMap> mLocalToFlatOverlapIndex; // (locali, localj) -> flatj

  // Surface connectivity: the surfaces with which each point overlaps
  std::vector<std::vector<Vector>> mSurfaceNormal; // (locali, flats) -> normal
  std::vector<NormalMap> mSurfaceFlatIndex; // (locali, normalasarray) -> flats
  std::vector<std::vector<int>> mVoidSurfaces; // (locali, flats) -> faceindexfori
  
  // Boundary information
  std::vector<int> mConstantBoundaryNodes; // list of the locali that are constantBoundaryNode
  std::vector<bool> mIsConstantBoundaryNode; // locali -> isConstantBoundaryNode
  std::vector<int> mNumConstantBoundaryNeighbors; // locali -> num neighbors that are constantBoundaryNode
  std::vector<int> mNumConstantBoundaryOverlapNeighbors; // locali -> num overlap neighbors that are constantBoundaryNode
  
  // All normals within this tolerance will be combined for indexing and integrals
  static constexpr double mRoundValue = 1.e8;
  static constexpr double mRoundValueInv = 1.e-8;
  
  // Scratch
  mutable ArrayDim mScratchArray;
};

} // end namespace Spheral

#include "FlatConnectivityInline.hh"

#endif
