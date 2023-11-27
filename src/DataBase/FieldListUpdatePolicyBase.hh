//---------------------------------Spheral++----------------------------------//
// FieldListUpdatePolicyBase -- Base/interface class for the policies on how 
// FieldList state variables are to be updated.
//
// Created by JMO, Sun Feb 13 20:50:53 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldListUpdatePolicyBase_hh__
#define __Spheral_FieldListUpdatePolicyBase_hh__

#include "UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldListUpdatePolicyBase: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  FieldListUpdatePolicyBase();
  explicit FieldListUpdatePolicyBase(const std::string& depend0);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1,
                            const std::string& depend2);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1, 
                            const std::string& depend2, 
                            const std::string& depend3);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1,
                            const std::string& depend2,
                            const std::string& depend3, 
                            const std::string& depend4);
  FieldListUpdatePolicyBase(const std::string& depend0, 
                            const std::string& depend1,
                            const std::string& depend2,
                            const std::string& depend3,
                            const std::string& depend4, 
                            const std::string& depend5);
  FieldListUpdatePolicyBase(const std::string& depend0,
                            const std::string& depend1,
                            const std::string& depend2, 
                            const std::string& depend3, 
                            const std::string& depend4,
                            const std::string& depend5,
                            const std::string& depend6);
  virtual ~FieldListUpdatePolicyBase() {};
  
private:
  //--------------------------- Private Interface ---------------------------//
  FieldListUpdatePolicyBase(const FieldListUpdatePolicyBase& rhs);
  FieldListUpdatePolicyBase& operator=(const FieldListUpdatePolicyBase& rhs);
};

}

#ifndef __GCCXML__
#include "FieldListUpdatePolicyBaseInline.hh"
#endif

#endif
