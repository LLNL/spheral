//---------------------------------Spheral++----------------------------------//
// FieldListView -- A list container for Fields.  This is how Spheral++ defines
//              "global" fields, or fields that extend over more than one
//              NodeList.
// A FieldListView can either just hold pointers to externally stored Fields, or
// copy the Fields to internal storage.
//
// Created by JMO, Sat Feb  5 12:57:58 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldSpace__FieldListView_hh__
#define __Spheral__FieldSpace__FieldListView_hh__

#include "Utilities/OpenMP_wrapper.hh"
#include "SphArray.hh"
#include "config.hh"
#include "FieldView.hh"

#include <vector>
#include <list>
#include <map>
#include <memory>

namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension, typename DataType>
class FieldListView :
    public chai::CHAICopyable{
public:
  //--------------------------- Public Interface ---------------------------//
  
  using FieldListType = FieldList<Dimension, DataType>;

  typedef typename Dimension::Scalar Scalar;

  //typedef Field<Dimension, DataType>* ElementType;
  //typedef Field<Dimension, DataType>* value_type;    // STL compatibility
  using ElementType = FieldView<Dimension, DataType>;
  using value_type  = ElementType;
  //using StorageType = std::vector<ElementType>;
  //using StorageType = ManagedVector<ElementType>;
  using StorageType = MVSmartRef<ElementType>;

  typedef typename StorageType::MV::iterator iterator;
  typedef typename StorageType::MV::const_iterator const_iterator;

  // Constructors.
  SPHERAL_HOST_DEVICE FieldListView();
  SPHERAL_HOST FieldListView(FieldListType const& rhs);
  SPHERAL_HOST_DEVICE FieldListView(FieldListView const& rhs);

  // Destructor.
  SPHERAL_HOST_DEVICE ~FieldListView();

  // Assignment.
  FieldListView& operator=(const FieldListView& rhs);
  FieldListView& operator=(const DataType& rhs);

  // Provide the standard iterators over the Fields.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // Index operator.
  ElementType operator[](const unsigned index);
  ElementType operator[](const unsigned index) const;

  ElementType at(const unsigned index);
  ElementType at(const unsigned index) const;

  // Provide a more primitive access to Field elements, based on the index of the Field
  // and the node index within that field.
  DataType& operator()(const unsigned int fieldIndex,
                       const unsigned int nodeIndex);
  DataType& operator()(const unsigned int fieldIndex,
                             const unsigned int nodeIndex) const;

  // Reproduce the standard Field operators for FieldListViews.
  void Zero();
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const double dataMin);
  void applyScalarMax(const double dataMax);

  //FieldListView operator+(const FieldListView& rhs) const;
  //FieldListView operator-(const FieldListView& rhs) const;

  //FieldListView& operator+=(const FieldListView& rhs);
  //FieldListView& operator-=(const FieldListView& rhs);

  //FieldListView operator+(const DataType& rhs) const;
  //FieldListView operator-(const DataType& rhs) const;

  //FieldListView& operator+=(const DataType& rhs);
  //FieldListView& operator-=(const DataType& rhs);

  //FieldListView<Dimension, DataType>& operator*=(const FieldListView<Dimension, Scalar>& rhs);
  //FieldListView<Dimension, DataType>& operator*=(const Scalar& rhs);

  //FieldListView<Dimension, DataType> operator/(const FieldListView<Dimension, Scalar>& rhs) const;
  //FieldListView<Dimension, DataType> operator/(const Scalar& rhs) const;

  //FieldListView<Dimension, DataType>& operator/=(const FieldListView<Dimension, Scalar>& rhs);
  //FieldListView<Dimension, DataType>& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const;
  DataType localMin() const;
  DataType localMax() const;

  //// Comparison operators (Field-Field element wise).
  bool operator==(const FieldListView& rhs) const;
  bool operator!=(const FieldListView& rhs) const;
  bool operator>(const FieldListView& rhs) const;
  bool operator<(const FieldListView& rhs) const;
  bool operator>=(const FieldListView& rhs) const;
  bool operator<=(const FieldListView& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator>(const DataType& rhs) const;
  bool operator<(const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  //// The number of fields in the FieldListView.
  unsigned numFields() const;
  unsigned size() const;

  ////----------------------------------------------------------------------------
  //// Methods to facilitate threaded computing
  //// Make a local thread copy of all the Fields
  //FieldListView<Dimension, DataType> threadCopy(const ThreadReduction reductionType = ThreadReduction::SUM,
  //                                          const bool copy = false);

  //// Same thing, with a "stack" object to simplify final reduction
  //FieldListView<Dimension, DataType> threadCopy(typename SpheralThreads<Dimension>::FieldListViewStack& stack,
  //                                          const ThreadReduction reductionType = ThreadReduction::SUM,
  //                                          const bool copy = false);

  //// Reduce the values in the FieldListView with the passed thread-local values.
  //void threadReduce() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  StorageType mFieldPtrs;

  //StorageType & mFieldViews() { return mFieldPtrs; }
  //StorageType const& mFieldViews() const { return mFieldPtrs; }
  typename StorageType::MV & mFieldViews() { return *(mFieldPtrs.get()); }
  typename StorageType::MV const& mFieldViews() const { return *(mFieldPtrs.get()); }

public:
  //// The master FieldListView if this is a thread copy.
  //FieldListView<Dimension, DataType>* threadMasterPtr;

};

} // namespace Spheral

#include "FieldListViewInline.hh"

#else

namespace Spheral {
  template<typename Dimension, typename DataType> class FieldListView;
}

#endif
