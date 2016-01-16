//---------------------------------Spheral++----------------------------------//
// DataBase -- The central point to store NodeLists and FieldLists, to serve the 
//             Spheral++ physics modules.
//
// Created by JMO, Sun Feb  6 13:44:49 PST 2000
//----------------------------------------------------------------------------//
#ifndef DataBase_HH
#define DataBase_HH

#ifndef __GCCXML__
#include <vector>
#include "boost/shared_ptr.hpp"
#else
#include "fakestl.hh"
#endif

#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Strength/SolidNodeList.hh"
#include "Field/NodeIterators.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace Kernel {
    template<typename Dimension> class TableKernel;
  }
}


namespace Spheral {
namespace DataBaseSpace {

template<typename Dimension>
class DataBase {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename std::vector<NodeSpace::NodeList<Dimension>*>::iterator NodeListIterator;
  typedef typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator ConstNodeListIterator;

  typedef typename std::vector<NodeSpace::FluidNodeList<Dimension>*>::iterator FluidNodeListIterator;
  typedef typename std::vector<NodeSpace::FluidNodeList<Dimension>*>::const_iterator ConstFluidNodeListIterator;

  typedef typename std::vector<SolidMaterial::SolidNodeList<Dimension>*>::iterator SolidNodeListIterator;
  typedef typename std::vector<SolidMaterial::SolidNodeList<Dimension>*>::const_iterator ConstSolidNodeListIterator;

  typedef NeighborSpace::ConnectivityMap<Dimension> ConnectivityMapType;
  typedef boost::shared_ptr<ConnectivityMapType> ConnectivityMapPtr;
  
  // It is convenient to be able to query the DataBase for the problem
  // dimensionality for Python.
  static int nDim;

  // Constructors.
  DataBase();

  // Destructor.
  ~DataBase();

  // Assignment.
  DataBase& operator=(const DataBase& rhs);

  // Number of NodeLists we have in the DataBase.
  int numNodeLists() const;
  int numFluidNodeLists() const;
  int numSolidNodeLists() const;

  // Numbers of nodes.
  int numInternalNodes() const;
  int numGhostNodes() const;
  int numNodes() const;

  int globalNumInternalNodes() const;
  int globalNumGhostNodes() const;
  int globalNumNodes() const;

  // Provide normal iterator methods over the DataBase NodeLists.
  NodeListIterator nodeListBegin();
  NodeListIterator nodeListEnd();

  ConstNodeListIterator nodeListBegin() const;
  ConstNodeListIterator nodeListEnd() const;

  FluidNodeListIterator fluidNodeListBegin();
  FluidNodeListIterator fluidNodeListEnd();

  ConstFluidNodeListIterator fluidNodeListBegin() const;
  ConstFluidNodeListIterator fluidNodeListEnd() const;

  NodeListIterator fluidNodeListAsNodeListBegin();
  NodeListIterator fluidNodeListAsNodeListEnd();

  ConstNodeListIterator fluidNodeListAsNodeListBegin() const;
  ConstNodeListIterator fluidNodeListAsNodeListEnd() const;

  SolidNodeListIterator solidNodeListBegin();
  SolidNodeListIterator solidNodeListEnd();

  ConstSolidNodeListIterator solidNodeListBegin() const;
  ConstSolidNodeListIterator solidNodeListEnd() const;

  NodeListIterator solidNodeListAsNodeListBegin();
  NodeListIterator solidNodeListAsNodeListEnd();

  ConstNodeListIterator solidNodeListAsNodeListBegin() const;
  ConstNodeListIterator solidNodeListAsNodeListEnd() const;

  // Provide NodeIterators to go over the elements of NodeLists/FieldLists.
  AllNodeIterator<Dimension> nodeBegin() const;
  AllNodeIterator<Dimension> nodeEnd() const;

  InternalNodeIterator<Dimension> internalNodeBegin() const;
  InternalNodeIterator<Dimension> internalNodeEnd() const;
  
  GhostNodeIterator<Dimension> ghostNodeBegin() const;
  GhostNodeIterator<Dimension> ghostNodeEnd() const;
  
  MasterNodeIterator<Dimension> masterNodeBegin() const;
  MasterNodeIterator<Dimension> masterNodeEnd() const;
  
  CoarseNodeIterator<Dimension> coarseNodeBegin() const;
  CoarseNodeIterator<Dimension> coarseNodeEnd() const;
  
  RefineNodeIterator<Dimension> refineNodeBegin() const;
  RefineNodeIterator<Dimension> refineNodeEnd() const;

  // Same iterator methods, but over FluidNodeLists.
  AllNodeIterator<Dimension> fluidNodeBegin() const;
  AllNodeIterator<Dimension> fluidNodeEnd() const;
  
  InternalNodeIterator<Dimension> fluidInternalNodeBegin() const;
  InternalNodeIterator<Dimension> fluidInternalNodeEnd() const;
  
  GhostNodeIterator<Dimension> fluidGhostNodeBegin() const;
  GhostNodeIterator<Dimension> fluidGhostNodeEnd() const;
  
  MasterNodeIterator<Dimension> fluidMasterNodeBegin() const;
  MasterNodeIterator<Dimension> fluidMasterNodeEnd() const;
  
  CoarseNodeIterator<Dimension> fluidCoarseNodeBegin() const;
  CoarseNodeIterator<Dimension> fluidCoarseNodeEnd() const;
  
  RefineNodeIterator<Dimension> fluidRefineNodeBegin() const;
  RefineNodeIterator<Dimension> fluidRefineNodeEnd() const;

  // Update the internal connectivity map.
  void updateConnectivityMap(const bool computeGhostConnectivity) const;
  void patchConnectivityMap(const FieldSpace::FieldList<Dimension, int>& flags,
                            const FieldSpace::FieldList<Dimension, int>& old2new) const;

  // Get the connectivity map.
  const ConnectivityMapType& connectivityMap() const;
  const ConnectivityMapType& connectivityMap(const bool computeGhostConnectivity) const;
  ConnectivityMapPtr connectivityMapPtr(const bool computeGhostConnectivity) const;

  // Methods to add, remove, and verify NodeLists.
  void appendNodeList(SolidMaterial::SolidNodeList<Dimension>& nodeList);
  void appendNodeList(NodeSpace::FluidNodeList<Dimension>& nodeList);
  void appendNodeList(NodeSpace::NodeList<Dimension>& nodeList);

  void deleteNodeList(SolidMaterial::SolidNodeList<Dimension>& nodeList);
  void deleteNodeList(NodeSpace::FluidNodeList<Dimension>& nodeList);
  void deleteNodeList(NodeSpace::NodeList<Dimension>& nodeList);

  bool haveNodeList(const NodeSpace::NodeList<Dimension>& nodeList) const;

  // Allow const access to the list of NodeList pointers.
  const std::vector<NodeSpace::NodeList<Dimension>*>& nodeListPtrs() const;
  const std::vector<NodeSpace::FluidNodeList<Dimension>*>& fluidNodeListPtrs() const;
  const std::vector<SolidMaterial::SolidNodeList<Dimension>*>& solidNodeListPtrs() const;

  // Provide convenience functions for manipulating the neighbor information
  // of the NodeLists.
  void setMasterNodeLists(const Vector& position,
                          const SymTensor& H) const;
  void setMasterFluidNodeLists(const Vector& position,
                               const SymTensor& H) const;

  void setRefineNodeLists(const Vector& position,
                          const SymTensor& H) const;
  void setRefineFluidNodeLists(const Vector& position,
                               const SymTensor& H) const;

  // Query methods which return "global" fields (FieldLists) for quantities
  // defined over all NodeLists.  These methods all build up FieldLists
  // referencing internally cached data in NodeLists, and so are fast.
  FieldSpace::FieldList<Dimension, Scalar> globalMass() const;
  FieldSpace::FieldList<Dimension, Vector> globalPosition() const;
  FieldSpace::FieldList<Dimension, Vector> globalVelocity() const;
  FieldSpace::FieldList<Dimension, SymTensor> globalHfield() const;
  FieldSpace::FieldList<Dimension, Scalar> globalWork() const;

  FieldSpace::FieldList<Dimension, Scalar> fluidMass() const;
  FieldSpace::FieldList<Dimension, Vector> fluidPosition() const;
  FieldSpace::FieldList<Dimension, Vector> fluidVelocity() const;
  FieldSpace::FieldList<Dimension, Scalar> fluidMassDensity() const;
  FieldSpace::FieldList<Dimension, Scalar> fluidSpecificThermalEnergy() const;
  FieldSpace::FieldList<Dimension, SymTensor> fluidHfield() const;
  FieldSpace::FieldList<Dimension, Scalar> fluidWork() const;

  // We can also return the node extent Fields stored in the Neighbor objects.
  FieldSpace::FieldList<Dimension, Vector> globalNodeExtent() const;
  FieldSpace::FieldList<Dimension, Vector> fluidNodeExtent() const;

  // These functions return FieldLists with Fields that have to be calculated and
  // stored, so they are more expensive.
  // These FieldLists represent state that must be computed on the fly.
  void globalHinverse(FieldSpace::FieldList<Dimension, SymTensor>& result) const;
  void fluidHinverse(FieldSpace::FieldList<Dimension, SymTensor>& result) const;
  void fluidPressure(FieldSpace::FieldList<Dimension, Scalar>& result) const;
  void fluidTemperature(FieldSpace::FieldList<Dimension, Scalar>& result) const;
  void fluidSoundSpeed(FieldSpace::FieldList<Dimension, Scalar>& result) const;
  void fluidVolume(FieldSpace::FieldList<Dimension, Scalar>& result) const;
  void fluidGamma(FieldSpace::FieldList<Dimension, Scalar>& result) const;
  void fluidLinearMomentum(FieldSpace::FieldList<Dimension, Vector>& result) const;
  void fluidTotalEnergy(FieldSpace::FieldList<Dimension, Scalar>& result) const;

  // Collect the number of neighbors for each node from the ConnectivityMap.
  FieldSpace::FieldList<Dimension, int> numNeighbors() const;

  // Create new FieldLists of size the number of NodeLists or FluidNodeLists.
  template<typename DataType>
  FieldSpace::FieldList<Dimension, DataType> newGlobalFieldList(const DataType value,
                                                                const typename FieldSpace::Field<Dimension, DataType>::FieldName name = "Unnamed Field") const;
  template<typename DataType>
  FieldSpace::FieldList<Dimension, DataType> newFluidFieldList(const DataType value,
                                                               const typename FieldSpace::Field<Dimension, DataType>::FieldName name = "Unnamed Field") const;
  template<typename DataType>
  FieldSpace::FieldList<Dimension, DataType> newSolidFieldList(const DataType value,
                                                               const typename FieldSpace::Field<Dimension, DataType>::FieldName name = "Unnamed Field") const;

  // Resize a FieldList to the number of NodeLists or FluidNodeLists.
  // Optionally we can also set all elements in the FieldList to the specified value.
  // Note that if the FieldList is resized it is reconstructed from scratch, so all elements
  // will get the new value regardless of resetValues.
  template<typename DataType>
  void resizeGlobalFieldList(FieldSpace::FieldList<Dimension, DataType>& fieldList,
                             const DataType value,
                             const typename FieldSpace::Field<Dimension, DataType>::FieldName name = "Unnamed Field",
                             const bool resetValues = true) const;
  template<typename DataType>
  void resizeFluidFieldList(FieldSpace::FieldList<Dimension, DataType>& fieldList,
                            const DataType value,
                            const typename FieldSpace::Field<Dimension, DataType>::FieldName name = "Unnamed Field",
                            const bool resetValues = true) const;
  template<typename DataType>
  void resizeSolidFieldList(FieldSpace::FieldList<Dimension, DataType>& fieldList,
                            const DataType value,
                            const typename FieldSpace::Field<Dimension, DataType>::FieldName name = "Unnamed Field",
                            const bool resetValues = true) const;

  // Get the maximum kernel extent for all NodeLists.
  Scalar maxKernelExtent() const;

  // Compute coordinates bounding all nodes in the DataBase.
  void boundingBox(Vector& xmin, Vector& xmax,
                   const bool ghost = true) const;
  void boundingBox(Vector& xmin, Vector& xmax,
                   const FieldSpace::FieldList<Dimension, int>& mask,
                   const bool ghost = true) const;

  // Return the local and global max sampling extents.
  void localSamplingBoundingVolume(Vector& centroid, double& radiusNodes, double& radiusSample,
                                   Vector& xminNodes, Vector& xmaxNodes,
                                   Vector& xminSample, Vector& xmaxSample) const;
  void globalSamplingBoundingVolume(Vector& centroid, double& radiusNodes, double& radiusSample,
                                    Vector& xminNodes, Vector& xmaxNodes,
                                    Vector& xminSample, Vector& xmaxSample) const;

  // Return the local min and max sampling extents for groupings of connected nodes.
  void localSamplingBoundingBoxes(std::vector<Vector>& xminima,
                                  std::vector<Vector>& xmaxima) const;
  void globalSamplingBoundingBoxes(std::vector<Vector>& xminima,
                                   std::vector<Vector>& xmaxima) const;

  // Make each fluid node list update their mass densities.
  void updateFluidMassDensity() const;

  // Provide a method to determine if the DataBase is in a minimally defined
  // valid state.
  bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  std::vector<NodeSpace::NodeList<Dimension>*> mNodeListPtrs;

  std::vector<NodeSpace::FluidNodeList<Dimension>*> mFluidNodeListPtrs;
  std::vector<NodeSpace::NodeList<Dimension>*> mFluidNodeListAsNodeListPtrs;

  std::vector<SolidMaterial::SolidNodeList<Dimension>*> mSolidNodeListPtrs;
  std::vector<NodeSpace::NodeList<Dimension>*> mSolidNodeListAsNodeListPtrs;

  mutable ConnectivityMapPtr mConnectivityMapPtr;
#endif

  // Prevent copying.
  DataBase(const DataBase& rhs);
};

}
}

#ifndef __GCCXML__
#include "DataBaseInline.hh"
#endif

#else

namespace Spheral {
  namespace DataBaseSpace {
    // Forward declaration.
    template<typename Dimension> class DataBase;
  }
}

#endif
