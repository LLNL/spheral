//---------------------------------Spheral++----------------------------------//
// DistributedBoundary -- Base class for distributed parallel boundary
// conditions, connecting NodeLists across parallel domains.
//
// Created by JMO, Thu Aug 23 21:34:32 PDT 2001
//----------------------------------------------------------------------------//

#ifndef DistributedBoundary_HH
#define DistributedBoundary_HH

#ifndef __GCCXML__
#include <vector>
#include <map>
#include <list>
#else
#include "fakestl.hh"
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "Boundary/Boundary.hh"

// Forward declarations.
namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace FieldSpace {
    template<typename Dimension> class FieldBase;
    template<typename Dimension, typename DataType> class Field;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class DistributedBoundary: public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  struct DomainBoundaryNodes {
    std::vector<int> sendNodes;
    std::vector<int> receiveNodes;
  };

  typedef std::map<int, DomainBoundaryNodes> DomainBoundaryNodeMap;
  typedef std::map<NodeSpace::NodeList<Dimension>*, DomainBoundaryNodeMap> NodeListDomainBoundaryNodeMap;

  // Constructors and destructors.
  DistributedBoundary();
  virtual ~DistributedBoundary();

  // Get this domain ID.
  int domainID() const;

  // Total number of domains.
  int numDomains() const;

  // Test if the given NodeList is communicated on this domain or not.
  bool communicatedNodeList(const NodeSpace::NodeList<Dimension>& nodeList) const;

  // Test if the given NodeList is communicated with the given domain.
  bool nodeListSharedWithDomain(const NodeSpace::NodeList<Dimension>& nodeList,
                                int neighborDomainID) const;

  // Allow read access to the communication information.
  const NodeListDomainBoundaryNodeMap& nodeListDomainBoundaryNodeMap() const;
  const DomainBoundaryNodeMap& domainBoundaryNodeMap(const NodeSpace::NodeList<Dimension>& nodeList) const;
  const DomainBoundaryNodes& domainBoundaryNodes(const NodeSpace::NodeList<Dimension>&,
                                                 int neighborDomainID) const;

  // Extract the current set of processors we're communicating with.
  void communicatedProcs(std::vector<int>& sendProcs,
			 std::vector<int>& recvProcs) const;

//   // Generic exchange method for Fields.
//   template<typename DataType>
//   void exchangeField(FieldSpace::Field<Dimension, DataType>& field) const;

//   // We also have a specialized version for Field<vector<double> >.
//   void exchangeField(FieldSpace::Field<Dimension, std::vector<double> >& field) const;

  // Non-blocking exchanges.
  template<typename DataType>
  void beginExchangeField(FieldSpace::Field<Dimension, DataType>& field) const;

  template<typename DataType>
  void beginExchangeFieldVariableSize(FieldSpace::Field<Dimension, DataType>& field) const;

  //**********************************************************************
  // Descendent Distributed Neighbors are required to provide the 
  // setGhostNodes method for DataBases.
  virtual void setAllGhostNodes(DataBaseSpace::DataBase<Dimension>& dataBase) = 0;

  // Override the Boundary method for culling ghost nodes.
  virtual void cullGhostNodes(const FieldSpace::FieldList<Dimension, int>& flagSet,
                              FieldSpace::FieldList<Dimension, int>& old2newIndexMap,
                              std::vector<int>& numNodesRemoved) override;

  // Use the given NodeList's neighbor object to select the ghost nodes.
  // This method should never be called for the DistributedBoundary!
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList) override;

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Scalar>>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Vector>>& field) const override;

  // Distributed boundaries don't have "violate" nodes, so override these
  // methods to no-ops.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList) override;
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList) override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const override;
  //**********************************************************************

  // Override the base method to finalize ghost boundaries.
  virtual void finalizeGhostBoundary() const override;

  // We do not want to use the parallel ghost nodes as generators.
  virtual bool meshGhostNodes() const override;

  // Unpack a packed set of Field values back into the Field.
  template<typename DataType>
  void unpackField(FieldSpace::Field<Dimension, DataType>& field,
                   const std::list< std::vector<char> >& packedValues) const;

  // Force the exchanges which have been registered to execute.
  void finalizeExchanges();

  // Update the control and ghost nodes of the base class.
  void setControlAndGhostNodes();

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Descendent classes get read/write access to the communication maps.
  NodeListDomainBoundaryNodeMap& accessNodeListDomainBoundaryNodeMap();
  DomainBoundaryNodeMap& accessDomainBoundaryNodeMap(const NodeSpace::NodeList<Dimension>& nodeList);
  DomainBoundaryNodes& accessDomainBoundaryNodes(const NodeSpace::NodeList<Dimension>&,
                                                 int neighborDomainID);

  // Convenience method to return an iterator to the DomainBoundaryNodes for the
  // given NodeList and domain pair.  If there isn't an entry for this pair already,
  // this method creates one.
  DomainBoundaryNodes& openDomainBoundaryNodes(NodeSpace::NodeList<Dimension>* nodeListPtr,
                                               const int domainID);

  // Inverse of above -- remove the DomainBoundaryNodes for a NodeList/procID pair.
  void removeDomainBoundaryNodes(NodeSpace::NodeList<Dimension>* nodeListPtr,
                                 const int domainID);

  // Override the Boundary method for clearing the maps.
  virtual void reset(const DataBaseSpace::DataBase<Dimension>& dataBase) override;

  // This handy helper method will build receive and ghost nodes on all each
  // domain based on send nodes that have already been filled in.
  void buildReceiveAndGhostNodes(const DataBaseSpace::DataBase<Dimension>& dataBase);

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  int mDomainID;
  NodeListDomainBoundaryNodeMap mNodeListDomainBoundaryNodeMap;
  
  // List of the fields that are currently backlogged for exchanging.
  mutable std::vector<FieldSpace::Field<Dimension, int>*> mIntExchangeFields;
  mutable std::vector<FieldSpace::Field<Dimension, Scalar>*> mScalarExchangeFields;
  mutable std::vector<FieldSpace::Field<Dimension, Vector>*> mVectorExchangeFields;
  mutable std::vector<FieldSpace::Field<Dimension, Tensor>*> mTensorExchangeFields;
  mutable std::vector<FieldSpace::Field<Dimension, SymTensor>*> mSymTensorExchangeFields;
  mutable std::vector<FieldSpace::Field<Dimension, ThirdRankTensor>*> mThirdRankTensorExchangeFields;
  mutable std::vector<FieldSpace::Field<Dimension, std::vector<Scalar>>*> mVectorScalarExchangeFields;
  mutable std::vector<FieldSpace::Field<Dimension, std::vector<Vector>>*> mVectorVectorExchangeFields;

  // Internal tag for MPI communiators.
  mutable int mMPIFieldTag;

#ifdef USE_MPI
  // Send/receive requests.
  mutable std::vector<MPI_Request> mSendRequests;
  mutable std::vector<MPI_Request> mRecvRequests;
#endif

#ifdef USE_MPI_DEADLOCK_DETECTION
  mutable std::vector<int> mSendProcIDs;
  mutable std::vector<int> mRecvProcIDs;
#endif

  // Buffers for holding send/receive data.
  typedef std::list< std::list< std::vector<char> > > CommBufferSet;
  mutable CommBufferSet mSendBuffers;
  mutable CommBufferSet mRecvBuffers;

  // Maps linking fields to their communication buffers.
  typedef std::map<const FieldSpace::FieldBase<Dimension>*, std::list< std::vector<char> >* > Field2BufferType;
  mutable Field2BufferType mField2SendBuffer;
  mutable Field2BufferType mField2RecvBuffer;
#endif

};

}
}

#ifndef __GCCXML__
#include "DistributedBoundaryInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    template<typename Dimension> class DistributedBoundary;
  }
}

#endif
