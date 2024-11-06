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

#include "DataOutput/registerWithRestart.hh"
#include "Utilities/ValueViewInterface.hh"

#include <string>
#include <list>
#include <vector>

namespace Spheral {

template<typename Dimension> class AllNodeIterator;
template<typename Dimension> class InternalNodeIterator;
template<typename Dimension> class GhostNodeIterator;
template<typename Dimension> class MasterNodeIterator;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class State;
template<typename Dimension> class Neighbor;
template<typename Dimension> class DataBase;
template<typename Dimension> class TableKernel;

template<typename Dimension> class FieldBase;
template<typename Dimension> class FieldBaseView;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

enum class NodeType {
  InternalNode = 0,
  GhostNode = 1
};

//VVI_IMPL_BEGIN
//template<typename Dimension> class FieldBase;
//template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class NodeList {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using FieldBaseIterator = typename std::vector<FieldBase<Dimension>*>::iterator;
  using const_FieldBaseIterator = typename std::vector<FieldBase<Dimension>*>::const_iterator;

private:
  //--------------------------- Private Interface ---------------------------//
  unsigned mNumNodes;
  unsigned mFirstGhostNode;

  std::string mName;

  // State fields.
  Field<Dimension, Scalar> mMass;
  Field<Dimension, Vector> mPositions;
  Field<Dimension, Vector> mVelocity;
  Field<Dimension, SymTensor> mH;

  // The work field is mutable.
  mutable Field<Dimension, Scalar> mWork;

  // Stuff for how H is handled.
  Scalar mhmin, mhmax, mhminratio;
  Scalar mNodesPerSmoothingScale;
  unsigned mMaxNumNeighbors;

  // List of fields that are defined over this NodeList.
  mutable std::vector<FieldBaseView<Dimension>> mFieldBaseList;
  Neighbor<Dimension>* mNeighborPtr;

  // Provide a dummy vector to buid NodeIterators against.
  std::vector<NodeList<Dimension>*> mDummyList;

  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  NodeList();
  NodeList(const NodeList& nodes);
  NodeList& operator=(const NodeList& rhs);

public:
  // Constructors
  explicit NodeList(std::string name,
                    const unsigned numInternal,
                    const unsigned numGhost,
                    const Scalar hmin = 1.0e-20,
                    const Scalar hmax = 1.0e20,
                    const Scalar hminratio = 0.1,
                    const Scalar nPerh = 2.01,
                    const unsigned maxNumNeighbors = 500);

  // Destructor
  virtual ~NodeList();

  // Access the name of the NodeList.
  std::string name() const;

  // Get or set the number of Nodes.
  unsigned numNodes() const;
  unsigned numInternalNodes() const;
  unsigned numGhostNodes() const;
  void numInternalNodes(unsigned size);
  void numGhostNodes(unsigned size);

  // Provide the standard NodeIterators over the nodes of this NodeList.
  AllNodeIterator<Dimension> nodeBegin() const;
  AllNodeIterator<Dimension> nodeEnd() const;
          
  InternalNodeIterator<Dimension> internalNodeBegin() const;
  InternalNodeIterator<Dimension> internalNodeEnd() const;
          
  GhostNodeIterator<Dimension> ghostNodeBegin() const;
  GhostNodeIterator<Dimension> ghostNodeEnd() const;
          
  MasterNodeIterator<Dimension> masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const;
  MasterNodeIterator<Dimension> masterNodeEnd() const;
          
  CoarseNodeIterator<Dimension> coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const;
  CoarseNodeIterator<Dimension> coarseNodeEnd() const;

  RefineNodeIterator<Dimension> refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const;
  RefineNodeIterator<Dimension> refineNodeEnd() const;

  // The NodeList state Fields.
  Field<Dimension, Scalar>& mass();
  Field<Dimension, Vector>& positions();
  Field<Dimension, Vector>& velocity();
  Field<Dimension, SymTensor>& Hfield();

  const Field<Dimension, Scalar>& mass() const;
  const Field<Dimension, Vector>& positions() const;
  const Field<Dimension, Vector>& velocity() const;
  const Field<Dimension, SymTensor>& Hfield() const;

  void mass(const Field<Dimension, Scalar>& m);
  void positions(const Field<Dimension, Vector>& r);
  void velocity(const Field<Dimension, Vector>& v);
  void Hfield(const Field<Dimension, SymTensor>& H);

  // Anyone can acces the work field as a mutable field.
  Field<Dimension, Scalar>& work() const;
  void work(const Field<Dimension, Scalar>& w);

  // These are quantities which are not stored, but can be computed.
  void Hinverse(Field<Dimension, SymTensor>& field) const;

  // Provide iterators over the set of FieldBases defined on this 
  // NodeList.
  FieldBaseIterator registeredFieldsBegin();
  FieldBaseIterator registeredFieldsEnd();

  const_FieldBaseIterator registeredFieldsBegin() const;
  const_FieldBaseIterator registeredFieldsEnd() const;

  // Provide methods to add and subtract Fields which are defined over a
  // NodeList.
  void registerField(FieldBaseView<Dimension> field) const;
  void unregisterField(FieldBaseView<Dimension> field) const;
  int numFields() const;
  bool haveField(const FieldBase<Dimension>& field) const;

  // NodeLists can contain ghost nodes (either communicated from neighbor
  // processors, or simply created for boundary conditions).
  NodeType nodeType(int i) const;
  unsigned firstGhostNode() const;

  // Access the neighbor object.
  Neighbor<Dimension>& neighbor() const;

  void registerNeighbor(Neighbor<Dimension>& neighbor);
  void unregisterNeighbor();

  // The target number of nodes per smoothing scale (for calculating the ideal H).
  Scalar nodesPerSmoothingScale() const;
  void nodesPerSmoothingScale(Scalar val);

  // The maximum number of neighbors we want to have (for calculating the ideal H).
  unsigned maxNumNeighbors() const;
  void maxNumNeighbors(unsigned val);

  // Allowed range of smoothing scales for use in calculating H.
  Scalar hmin() const;
  void hmin(Scalar val);

  Scalar hmax() const;
  void hmax(Scalar val);

  Scalar hminratio() const;
  void hminratio(Scalar val);

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
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Some operators.
  bool operator==(const NodeList& rhs) const;
  bool operator!=(const NodeList& rhs) const;

};

//VVI_IMPL_END
//
//#ifdef VVI_ENABLED
//
//template<typename Dimension>
//class NodeList;
//
//#define NodeListView__(code) PTR_VIEW_METACLASS_DECL((NodeList<Dimension>), (NodeListView), (vvimpl::NodeList<Dimension>), (code))
//#define NodeList__(code) PTR_VALUE_METACLASS_DECL((NodeList), (NodeListView<Dimension>), (UNPACK code))
//
//
//template<typename Dimension>
//class NodeListView__(
//public:
//  // NodeList inherits from NodeListBase so we would like to be able to implicitly
//  // upcast a fieldview object to a fieldbaseview.
//  //VVI_UPCAST_CONVERSION_OP((NodeListBaseView<Dimension>))
//);
//
//
//template<typename Dimension>
//class NodeList__((
//public:
//  using Scalar = typename ImplType::Scalar;
//  using Vector = typename ImplType::Vector;
//  using Tensor = typename ImplType::Tensor;
//  using SymTensor = typename ImplType::SymTensor;
//
//  using FieldBaseIterator = typename ImplType::FieldBaseIterator;
//  using const_FieldBaseIterator = typename ImplType::const_FieldBaseIterator;
//// Access the name of the NodeList.
//  std::string name() const { return VVI_IMPL_INST().name(); }
//
//  // Get or set the number of Nodes.
//  unsigned numNodes() const { return VVI_IMPL_INST().numNodes(); }
//  unsigned numInternalNodes() const { return VVI_IMPL_INST().numInternalNodes(); }
//  unsigned numGhostNodes() const { return VVI_IMPL_INST().numGhostNodes(); }
//  void numInternalNodes(unsigned size) { VVI_IMPL_INST().numInternalNodes(size); }
//  void numGhostNodes(unsigned size) { VVI_IMPL_INST().numGhostNodes(size); }
//
//  // Provide the standard NodeIterators over the nodes of this NodeList.
//  AllNodeIterator<Dimension> nodeBegin() const { return VVI_IMPL_INST().nodeBegin(); }
//  AllNodeIterator<Dimension> nodeEnd() const { return VVI_IMPL_INST().nodeEnd(); }
//          
//  InternalNodeIterator<Dimension> internalNodeBegin() const { return VVI_IMPL_INST().internalNodeBegin(); }
//  InternalNodeIterator<Dimension> internalNodeEnd() const { return VVI_IMPL_INST().internalNodeEnd(); }
//          
//  GhostNodeIterator<Dimension> ghostNodeBegin() const { return VVI_IMPL_INST().ghostNodeBegin(); }
//  GhostNodeIterator<Dimension> ghostNodeEnd() const { return VVI_IMPL_INST().ghostNodeEnd(); }
//          
//  MasterNodeIterator<Dimension> masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const { return VVI_IMPL_INST().masterNodeBegin(masterLists); }
//  MasterNodeIterator<Dimension> masterNodeEnd() const { return VVI_IMPL_INST().masterNodeEnd(); }
//          
//  CoarseNodeIterator<Dimension> coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const { return VVI_IMPL_INST().coarseNodeBegin(coarseNeighbors); }
//  CoarseNodeIterator<Dimension> coarseNodeEnd() const { return VVI_IMPL_INST().coarseNodeEnd(); }
//
//  RefineNodeIterator<Dimension> refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const { return VVI_IMPL_INST().refineNodeBegin(refineNeighbors); }
//  RefineNodeIterator<Dimension> refineNodeEnd() const { return VVI_IMPL_INST().refineNodeEnd(); }
//
//  // The NodeList state Fields.
//  Field<Dimension, Scalar>& mass() { return VVI_IMPL_INST().mass(); }
//  Field<Dimension, Vector>& positions() { return VVI_IMPL_INST().positions(); }
//  Field<Dimension, Vector>& velocity() { return VVI_IMPL_INST().velocity(); }
//  Field<Dimension, SymTensor>& Hfield() { return VVI_IMPL_INST().Hfield(); }
//
//  const Field<Dimension, Scalar>& mass() const { return VVI_IMPL_INST().mass(); }
//  const Field<Dimension, Vector>& positions() const { return VVI_IMPL_INST().positions(); }
//  const Field<Dimension, Vector>& velocity() const { return VVI_IMPL_INST().velocity(); }
//  const Field<Dimension, SymTensor>& Hfield() const { return VVI_IMPL_INST().Hfield(); }
//
//  void mass(const Field<Dimension, Scalar>& m) { return VVI_IMPL_INST().mass(m); }
//  void positions(const Field<Dimension, Vector>& r) { VVI_IMPL_INST().positions(r); }
//  void velocity(const Field<Dimension, Vector>& v) { VVI_IMPL_INST().velocity(v); }
//  void Hfield(const Field<Dimension, SymTensor>& H) { VVI_IMPL_INST().Hfield(H); }
//
//  // Anyone can acces the work field as a mutable field.
//  Field<Dimension, Scalar>& work() const { return VVI_IMPL_INST().work(); }
//  void work(const Field<Dimension, Scalar>& w) { VVI_IMPL_INST().work(w); }
//
//  // These are quantities which are not stored, but can be computed.
//  void Hinverse(Field<Dimension, SymTensor>& field) const { VVI_IMPL_INST().Hinverse(field); }
//
//  // Provide iterators over the set of FieldBases defined on this 
//  // NodeList.
//  FieldBaseIterator registeredFieldsBegin() { return VVI_IMPL_INST().registeredFieldsBegin(); }
//  FieldBaseIterator registeredFieldsEnd() { return VVI_IMPL_INST().registeredFieldsEnd(); }
//
//  const_FieldBaseIterator registeredFieldsBegin() const { return VVI_IMPL_INST().registeredFieldsBegin(); }
//  const_FieldBaseIterator registeredFieldsEnd() const { return VVI_IMPL_INST().registeredFieldsEnd(); }
//
//  // Provide methods to add and subtract Fields which are defined over a
//  // NodeList.
//  void registerField(Spheral::vvimpl::FieldBase<Dimension>& field) const { VVI_IMPL_INST().registerField(field); }
//  void unregisterField(Spheral::vvimpl::FieldBase<Dimension>& field) const { VVI_IMPL_INST().unregisterField(field); }
//  int numFields() const { return VVI_IMPL_INST().numFields(); }
//  bool haveField(const FieldBase<Dimension>& field) const { return VVI_IMPL_INST().haveField(field); }
//
//  // NodeLists can contain ghost nodes (either communicated from neighbor
//  // processors, or simply created for boundary conditions).
//  NodeType nodeType(int i) const { return VVI_IMPL_INST().nodeType(i); }
//  unsigned firstGhostNode() const { return VVI_IMPL_INST().firstGhostNode(); }
//
//  // Access the neighbor object.
//  Neighbor<Dimension>& neighbor() const { return VVI_IMPL_INST().neighbor(); }
//
//  void registerNeighbor(Neighbor<Dimension>& neighbor) { VVI_IMPL_INST().registerNeighbor(neighbor); }
//  void unregisterNeighbor() { VVI_IMPL_INST().unregisterNeighbor(); }
//
//
////  // Custom Ctor, note we need to create the underlying implementation 
////  // object on ctor of value interfaces.
//  Scalar nodesPerSmoothingScale() const { return VVI_IMPL_INST().nodesPerSmoothingScale(); }
//  void nodesPerSmoothingScale(Scalar val) { VVI_IMPL_INST().nodesPerSmoothingScale(val); }
//
//  unsigned maxNumNeighbors() const { return VVI_IMPL_INST().maxNumNeighbors(); }
//  void maxNumNeighbors(unsigned val) { VVI_IMPL_INST().maxNumNeighbors(val); }
//
//  Scalar hmin() const { return VVI_IMPL_INST().hmin(); }
//  void hmin(Scalar val) { VVI_IMPL_INST().hmin(val); }
//
//  Scalar hmax() const { return VVI_IMPL_INST().hmax(); }
//  void hmax(Scalar val) { VVI_IMPL_INST().hmax(val); }
//
//  Scalar hminratio() const { return VVI_IMPL_INST().hminratio(); }
//  void hminratio(Scalar val) { VVI_IMPL_INST().hminratio(val); }
//
//  //****************************************************************************
//  // Methods for adding/removing individual nodes to/from the NodeList
//  // (including all Field information.  These methods are primarily useful
//  // for redistributing Nodes between parallel domains.
//  //void deleteNodes(const std::vector<int>& nodeIDs) { VVI_IMPL_INST().deleteNodes(nodeIDs); }
//  //std::list< std::vector<char> >  packNodeFieldValues(const std::vector<int>& nodeIDs) const { return VVI_IMPL_INST().packNodeFieldValues(nodeIDs); }
//  //void appendInternalNodes(const int numNewNodes, const std::list< std::vector<char> >& packedFieldValues) { VVI_IMPL_INST().appendInternalNodes(numNewNodes, packedFieldValues); }
//
//  //// A related method for reordering the nodes.
//  //void reorderNodes(const std::vector<int>& newOrdering) { VVI_IMPL_INST().reorderNodes(newOrdering); }
//  ////****************************************************************************
//
//  ////****************************************************************************
//  //// Methods required for restarting.
//  //// Dump and restore the NodeList state.
//  //std::string label() const { return "NodeList";  { return VVI_IMPL_INST().label(); }
//  //void dumpState(FileIO& file, const std::string& pathName) const { VVI_IMPL_INST().dumpState(file, pathName); }
//  //void restoreState(const FileIO& file, const std::string& pathName) { VVI_IMPL_INST().restoreState(file, pathName); }
//  ////****************************************************************************
//
//  //// Some operators.
//  //bool operator==(const NodeList& rhs) const { return VVI_IMPL_INST().operator==(rhs); }
//  //bool operator!=(const NodeList& rhs) const { return VVI_IMPL_INST().operator!=(rhs); }
//  ));
//  
//#endif // !defined(SPHERAL_ENABLE_VVI)


}

#include "NodeListInline.hh"

#endif
