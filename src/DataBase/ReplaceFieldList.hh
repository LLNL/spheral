//---------------------------------Spheral++----------------------------------//
// ReplaceFieldList -- An implementation of FieldListUpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceFieldList_hh__
#define __Spheral_ReplaceFieldList_hh__

#include "FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class ReplaceFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  ReplaceFieldList();
  explicit ReplaceFieldList(const std::string& depend0);
  ReplaceFieldList(const std::string& depend0, const std::string& depend1);
  ReplaceFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2);
  ReplaceFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3);
  ReplaceFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4);
  ReplaceFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5);
  virtual ~ReplaceFieldList();
  
  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // An alternate method to be called when you want to specify that the "Replace" information
  // in the derivatives is invalid, and instead the value should be treated as a time advancement
  // algorithm instead.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return "new "; }

private:
  //--------------------------- Private Interface ---------------------------//
  ReplaceFieldList(const ReplaceFieldList& rhs);
  ReplaceFieldList& operator=(const ReplaceFieldList& rhs);
};

}

#endif
