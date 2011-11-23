//---------------------------------Spheral++----------------------------------//
// NodeList -- An abstract base class for the NodeLists.
//
// We will define here the basic functionality we expect all NodeLists to 
// provide.
//
// Created by JMO, Tue Jun  8 23:58:36 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeList__
#define __Spheral_NodeList__

#include <string>

#ifndef __GCCXML__
#include <list>
#include <vector>
#include "DataOutput/registerWithRestart.hh"
#else
#include "fakestl.hh"
#endif

namespace Spheral {
  template<typename Dimension> class AllNodeIterator;
  template<typename Dimension> class InternalNodeIterator;
  template<typename Dimension> class GhostNodeIterator;
  template<typename Dimension> class MasterNodeIterator;
  template<typename Dimension> class CoarseNodeIterator;
  template<typename Dimension> class RefineNodeIterator;
  template<typename Dimension> class State;
  namespace NeighborSpace {
    template<typename Dimension> class Neighbor;
  }
  namespace FileIOSpace {
    class FileIO;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension> class FieldBase;
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
}

namespace Spheral {
namespace NodeSpace {

enum NodeType {
  InternalNode = 0,
  GhostNode = 1
};

template<typename Dimension>
class NodeList {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename std::vector<FieldSpace::FieldBase<Dimension>*>::iterator FieldBaseIterator;
  typedef typename std::vector<FieldSpace::FieldBase<Dimension>*>::const_iterator const_FieldBaseIterator;

  // Constructors
  explicit NodeList(std::string name,
                    const int numInternal,
                    const int numGhost,
                    const Scalar hmin = 1.0e-20,
                    const Scalar hmax = 1.0e20,
                    const Scalar hminratio = 0.1,
                    const Scalar nPerh = 2.01,
                    const int maxNumNeighbors = 500);

  // Destructor
  virtual ~NodeList();

  // Access the name of the NodeList.
  std::string name() const;

  // Get or set the number of Nodes.
  int numNodes() const;
  int numInternalNodes() const;
  int numGhostNodes() const;
  void numInternalNodes(int size);
  void numGhostNodes(int size);

  // Provide the standard NodeIterators over the nodes of this NodeList.
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

  // The NodeList state Fields.
  FieldSpace::Field<Dimension, Scalar>& mass();
  FieldSpace::Field<Dimension, Vector>& positions();
  FieldSpace::Field<Dimension, Vector>& velocity();
  FieldSpace::Field<Dimension, SymTensor>& Hfield();

  const FieldSpace::Field<Dimension, Scalar>& mass() const;
  const FieldSpace::Field<Dimension, Vector>& positions() const;
  const FieldSpace::Field<Dimension, Vector>& velocity() const;
  const FieldSpace::Field<Dimension, SymTensor>& Hfield() const;

  void mass(const FieldSpace::Field<Dimension, Scalar>& m);
  void positions(const FieldSpace::Field<Dimension, Vector>& r);
  void velocity(const FieldSpace::Field<Dimension, Vector>& v);
  void Hfield(const FieldSpace::Field<Dimension, SymTensor>& H);

  // Anyone can acces the work field as a mutable field.
  FieldSpace::Field<Dimension, Scalar>& work() const;
  void work(const FieldSpace::Field<Dimension, Scalar>& w);

  // These are quantities which are not stored, but can be computed.
  void Hinverse(FieldSpace::Field<Dimension, SymTensor>& field) const;

  // Provide iterators over the set of FieldBases defined on this 
  // NodeList.
  FieldBaseIterator registeredFieldsBegin();
  FieldBaseIterator registeredFieldsEnd();

  const_FieldBaseIterator registeredFieldsBegin() const;
  const_FieldBaseIterator registeredFieldsEnd() const;

  // Provide methods to add and subtract Fields which are defined over a
  // NodeList.
  void registerField(FieldSpace::FieldBase<Dimension>& field) const;
  void unregisterField(FieldSpace::FieldBase<Dimension>& field) const;
  int numFields() const;
  bool haveField(const FieldSpace::FieldBase<Dimension>& field) const;

  // NodeLists can contain ghost nodes (either communicated from neighbor
  // processors, or simply created for boundary conditions).
  NodeType nodeType(int i) const;
  int firstGhostNode() const;

  // Access the neighbor object.
  NeighborSpace::Neighbor<Dimension>& neighbor() const;

  void registerNeighbor(NeighborSpace::Neighbor<Dimension>& neighbor);
  void unregisterNeighbor();

  // The target number of nodes per smoothing scale (for calculating the ideal H).
  Scalar nodesPerSmoothingScale() const;
  void nodesPerSmoothingScale(const Scalar val);

  // The maximum number of neighbors we want to have (for calculating the ideal H).
  int maxNumNeighbors() const;
  void maxNumNeighbors(const int val);

  // Allowed range of smoothing scales for use in calculating H.
  Scalar hmin() const;
  void hmin(const Scalar val);

  Scalar hmax() const;
  void hmax(const Scalar val);

  Scalar hminratio() const;
  void hminratio(const Scalar val);

  //****************************************************************************
  // Methods for adding/removing individual nodes to/from the NodeList
  // (including all Field information.  These methods are primarily useful
  // for redistributing Nodes between parallel domains.
  virtual void deleteNodes(const std::vector<int>& nodeIDs);
  virtual std::list< std::vector<char> >  packNodeFieldValues(const std::vector<int>& nodeIDs) const;
  virtual void appendInternalNodes(const int numNewNodes,
                                   const std::list< std::vector<char> >& packedFieldValues);

  // A related method for reordering the nodes.
  virtual void reorderNodes(const std::vector<int>& newOrdering);
  //****************************************************************************

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "NodeList"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Some operators.
  bool operator==(const NodeList& rhs) const;
  bool operator!=(const NodeList& rhs) const;

protected:
  //--------------------------- Protected Interface ---------------------------//

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  int mNumNodes;
  int mFirstGhostNode;

  std::string mName;

  // State fields.
  FieldSpace::Field<Dimension, Scalar> mMass;
  FieldSpace::Field<Dimension, Vector> mPositions;
  FieldSpace::Field<Dimension, Vector> mVelocity;
  FieldSpace::Field<Dimension, SymTensor> mH;

  // The work field is mutable.
  mutable FieldSpace::Field<Dimension, Scalar> mWork;

  // Stuff for how H is handled.
  Scalar mhmin, mhmax, mhminratio;
  Scalar mNodesPerSmoothingScale;
  int mMaxNumNeighbors;

  // List of fields that are defined over this NodeList.
  mutable std::vector<FieldSpace::FieldBase<Dimension>*> mFieldBaseList;
  NeighborSpace::Neighbor<Dimension>* mNeighborPtr;

  // Provide a dummy vector to buid NodeIterators against.
  std::vector<NodeList<Dimension>*> mDummyList;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif

  // No default constructor, copying, or assignment.
  NodeList();
  NodeList(const NodeList& nodes);
  NodeList& operator=(const NodeList& rhs);
};

}
}

#ifndef __GCCXML__
#include "NodeListInline.hh"
#endif

#else
// Forward declare the NodeList class.
namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
}

#endif
