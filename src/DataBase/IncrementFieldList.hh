//---------------------------------Spheral++----------------------------------//
// IncrementFieldList -- An implementation of FieldListUpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementFieldList_hh__
#define __Spheral_IncrementFieldList_hh__

#include "FieldListUpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class IncrementFieldList: public FieldListUpdatePolicyBase<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename FieldListUpdatePolicyBase<Dimension, ValueType>::KeyType KeyType;

  // Constructors, destructor.
  IncrementFieldList(const bool wildCardDerivs = false);
  explicit IncrementFieldList(const std::string& depend0,
                              const bool wildCardDerivs = false);
  IncrementFieldList(const std::string& depend0, const std::string& depend1,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4,
                     const bool wildCardDerivs = false);
  IncrementFieldList(const std::string& depend0, const std::string& depend1, const std::string& depend2, const std::string& depend3, const std::string& depend4, const std::string& depend5,
                     const bool wildCardDerivs = false);
  virtual ~IncrementFieldList();
  
  // Overload the methods describing how to update FieldLists.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return "delta "; }

  // Flip whether we try to find multiple registered increment fields.
  bool wildCardDerivs() const;
  void wildCardDerivs(const bool val);

private:
  //--------------------------- Private Interface ---------------------------//
  IncrementFieldList(const IncrementFieldList& rhs);
  IncrementFieldList& operator=(const IncrementFieldList& rhs);

  // Flag for looking for multiple increment derivatives.
  bool mWildCardDerivs;
};

}

#include "IncrementFieldListInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension, typename ValueType> class IncrementFieldList;
}

#endif
