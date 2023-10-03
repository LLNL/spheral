//---------------------------------Spheral++----------------------------------//
// ReplaceFieldList -- An implementation of FieldListUpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceFieldList_hh__
#define __Spheral_ReplaceFieldList_hh__

#include "FieldListUpdatePolicyBase.hh"
#include "ReplaceState.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class ReplaceFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  explicit ReplaceFieldList(const FieldList<Dimension, ValueType>& fieldList);
  ReplaceFieldList(const FieldList<Dimension, ValueType>& fieldList,
                   const std::string& depend0);
  ReplaceFieldList(const FieldList<Dimension, ValueType>& fieldList,
                   const std::string& depend0, const std::string& depend1);
  ReplaceFieldList(const FieldList<Dimension, ValueType>& fieldList,
                   const std::string& depend0, const std::string& depend1, const std::string& depend2);
  ReplaceFieldList(const FieldList<Dimension, ValueType>& fieldList,
                   const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  ReplaceFieldList(const FieldList<Dimension, ValueType>& fieldList,
                   const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  ReplaceFieldList(const FieldList<Dimension, ValueType>& fieldList,
                   const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~ReplaceFieldList();
  
  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return ReplaceState<Dimension, ValueType>::prefix(); }

private:
  //--------------------------- Private Interface ---------------------------//
  ReplaceFieldList(const ReplaceFieldList& rhs);
  ReplaceFieldList& operator=(const ReplaceFieldList& rhs);
};

}

#include "ReplaceFieldListInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class ReplaceFieldList;
}

#endif
