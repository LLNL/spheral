//---------------------------------Spheral++----------------------------------//
// FieldList -- A list container for Fields.  This is how Spheral++ defines
//              "global" fields, or fields that extend over more than one
//              NodeList.
// A FieldList can either just hold pointers to externally stored Fields, or
// copy the Fields to internal storage.
//
// Created by JMO, Sat Feb  5 12:57:58 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldSpace__FieldList_hh__
#define __Spheral__FieldSpace__FieldList_hh__

#ifndef __GCCXML__
#include <vector>
#include <list>
#include <map>
#include "boost/shared_ptr.hpp"
#else
#include "fakestl.hh"
#endif

#include "FieldListBase.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class NodeIteratorBase;
  template<typename Dimension> class AllNodeIterator;
  template<typename Dimension> class InternalNodeIterator;
  template<typename Dimension> class GhostNodeIterator;
  template<typename Dimension> class MasterNodeIterator;
  template<typename Dimension> class CoarseNodeIterator;
  template<typename Dimension> class RefineNodeIterator;
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace FieldSpace {

// An enum for selecting how Fields are stored in FieldLists.
enum FieldStorageType {
  Reference = 0,
  Copy = 1
};

template<typename Dimension, typename DataType>
class FieldList: public FieldListBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  
  typedef DataType FieldDataType;

  typedef FieldBase<Dimension>* BaseElementType;
  typedef Field<Dimension, DataType>* ElementType;
  typedef Field<Dimension, DataType>* value_type;    // STL compatibility
  typedef std::vector<ElementType> StorageType;

  typedef typename StorageType::iterator iterator;
  typedef typename StorageType::const_iterator const_iterator;
  typedef typename StorageType::reverse_iterator reverse_iterator;
  typedef typename StorageType::const_reverse_iterator const_reverse_iterator;

  typedef std::vector<DataType> CacheElementsType;
  typedef typename CacheElementsType::iterator cache_iterator;
  typedef typename CacheElementsType::const_iterator const_cache_iterator;

  // Constructors.
  FieldList();
  explicit FieldList(FieldStorageType aStorageType);
  FieldList(const FieldList& rhs);

  // Destructor.
  ~FieldList();

  // Assignment.
  FieldList& operator=(const FieldList& rhs);
  FieldList& operator=(const DataType& rhs);

  // Access the storage type of the field list.
  FieldStorageType storageType() const;

  // Force the Field storage to be Copy.
  void copyFields();

  // Test if the given field (or NodeList) is part of a FieldList.
  bool haveField(const Field<Dimension, DataType>& field) const;
  bool haveNodeList(const NodeSpace::NodeList<Dimension>& nodeList) const;

  // Force the Field members of this FieldList to be equal to those of
  // another FieldList.
  void assignFields(const FieldList& fieldList);

  // Convenience methods to add and delete Fields.
  void appendField(const Field<Dimension, DataType>& field);
  void deleteField(const Field<Dimension, DataType>& field);

  // Construct a new field and add it to the FieldList.
  // Note this only makes sense when we're storing fields as copies!
  void appendNewField(const typename FieldSpace::Field<Dimension, DataType>::FieldName name,
                      const NodeSpace::NodeList<Dimension>& nodeList,
                      const DataType value);

  // Provide the standard iterators over the Fields.
  iterator begin();
  iterator end();
  reverse_iterator rbegin();
  reverse_iterator rend();

  const_iterator begin() const;
  const_iterator end() const;
  const_reverse_iterator rbegin() const;
  const_reverse_iterator rend() const;

  // Iterators over FieldBase* required by base class.
  virtual typename FieldListBase<Dimension>::iterator begin_base();
  virtual typename FieldListBase<Dimension>::iterator end_base();
  virtual typename FieldListBase<Dimension>::reverse_iterator rbegin_base();
  virtual typename FieldListBase<Dimension>::reverse_iterator rend_base();

  virtual typename FieldListBase<Dimension>::const_iterator begin_base() const;
  virtual typename FieldListBase<Dimension>::const_iterator end_base() const;
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rbegin_base() const;
  virtual typename FieldListBase<Dimension>::const_reverse_iterator rend_base() const;

  // Index operator.
  ElementType operator[](const unsigned index);
  ElementType operator[](const unsigned index) const;

  ElementType at(const unsigned index);
  ElementType at(const unsigned index) const;

  // Return an iterator to the Field associated with the given NodeList.
  iterator fieldForNodeList(const NodeSpace::NodeList<Dimension>& nodeList);
  const_iterator fieldForNodeList(const NodeSpace::NodeList<Dimension>& nodeList) const;

  // Provide access to the Field elements via NodeIterators.
  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  // Provide a more primitive access to Field elements, based on the index of the Field
  // and the node index within that field.
  DataType& operator()(const unsigned int fieldIndex,
                       const unsigned int nodeIndex);
  const DataType& operator()(const unsigned int fieldIndex,
                             const unsigned int nodeIndex) const;

  // Return the interpolated value of the FieldList at a position.
  DataType operator()(const Vector& position,
                      const KernelSpace::TableKernel<Dimension>& W) const;

  // Provide NodeIterators on the elements of the FieldList.
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

  // Provide a convenience function for setting the neighbor node information
  // for all the NodeList in this FieldList.
  void setMasterNodeLists(const Vector& r, const SymTensor& H) const;
  void setMasterNodeLists(const Vector& r) const;

  void setRefineNodeLists(const Vector& r, const SymTensor& H) const;
  void setRefineNodeLists(const Vector& r) const;

  // Reproduce the standard Field operators for FieldLists.
  void Zero();
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const double dataMin);
  void applyScalarMax(const double dataMax);

  FieldList operator+(const FieldList& rhs) const;
  FieldList operator-(const FieldList& rhs) const;

  FieldList& operator+=(const FieldList& rhs);
  FieldList& operator-=(const FieldList& rhs);

  FieldList operator+(const DataType& rhs) const;
  FieldList operator-(const DataType& rhs) const;

  FieldList& operator+=(const DataType& rhs);
  FieldList& operator-=(const DataType& rhs);

  FieldList<Dimension, DataType>& operator*=(const FieldList<Dimension, Scalar>& rhs);
  FieldList<Dimension, DataType>& operator*=(const Scalar& rhs);

  FieldList<Dimension, DataType> operator/(const FieldList<Dimension, Scalar>& rhs) const;
  FieldList<Dimension, DataType> operator/(const Scalar& rhs) const;

  FieldList<Dimension, DataType>& operator/=(const FieldList<Dimension, Scalar>& rhs);
  FieldList<Dimension, DataType>& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const;
  DataType localMin() const;
  DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  bool operator==(const FieldList& rhs) const;
  bool operator!=(const FieldList& rhs) const;
  bool operator>(const FieldList& rhs) const;
  bool operator<(const FieldList& rhs) const;
  bool operator>=(const FieldList& rhs) const;
  bool operator<=(const FieldList& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator>(const DataType& rhs) const;
  bool operator<(const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  // The number of fields in the FieldList.
  int numFields() const;
  int size() const;

  // The number of nodes in the FieldList.
  int numNodes() const;
  
  // The number of internal nodes in the FieldList.
  int numInternalNodes() const;
  
  // The number of ghost nodes in the FieldList.
  int numGhostNodes() const;

  // Get the NodeLists this FieldList is defined on.
  const std::vector<NodeSpace::NodeList<Dimension>*>& nodeListPtrs() const;

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  typedef std::list<boost::shared_ptr<Field<Dimension, DataType> > > FieldCacheType;
  typedef std::map<const NodeSpace::NodeList<Dimension>*, int> HashMapType;

  std::vector<ElementType> mFieldPtrs;
  std::vector<BaseElementType> mFieldBasePtrs;
  FieldCacheType mFieldCache;
  FieldStorageType mStorageType;

  // Maintain a vector of the NodeLists this FieldList is defined in order to
  // construct NodeIterators.
  std::vector<NodeSpace::NodeList<Dimension>*> mNodeListPtrs;
  HashMapType mNodeListIndexMap;

  // Internal method to build the NodeListIndexMap from scratch.
  void buildNodeListIndexMap();
#endif

};

}
}

#ifndef __GCCXML__
#include "FieldListInline.hh"
#endif

#else

namespace Spheral {
  namespace FieldSpace {
    // Forward declaration.
    template<typename Dimension, typename DataType> class FieldList;
  }
}

#endif
