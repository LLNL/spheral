//---------------------------------Spheral++----------------------------------//
// FieldUpdatePolicyBase -- Base/interface class for the policies on how 
// Field state variables are to be updated.
//
// Created by JMO, Sun Feb 13 20:50:53 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldUpdatePolicyBase_hh__
#define __Spheral_FieldUpdatePolicyBase_hh__

#include "UpdatePolicyBase.hh"

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldUpdatePolicyBase: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  FieldUpdatePolicyBase();
  explicit FieldUpdatePolicyBase(const std::string& depend0);
  FieldUpdatePolicyBase(const std::string& depend0, 
                        const std::string& depend1);
  FieldUpdatePolicyBase(const std::string& depend0, 
                        const std::string& depend1,
                        const std::string& depend2);
  FieldUpdatePolicyBase(const std::string& depend0, 
                        const std::string& depend1, 
                        const std::string& depend2, 
                        const std::string& depend3);
  FieldUpdatePolicyBase(const std::string& depend0, 
                        const std::string& depend1,
                        const std::string& depend2,
                        const std::string& depend3, 
                        const std::string& depend4);
  FieldUpdatePolicyBase(const std::string& depend0, 
                        const std::string& depend1,
                        const std::string& depend2,
                        const std::string& depend3,
                        const std::string& depend4, 
                        const std::string& depend5);
  FieldUpdatePolicyBase(const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2, 
                        const std::string& depend3, 
                        const std::string& depend4,
                        const std::string& depend5,
                        const std::string& depend6);
  virtual ~FieldUpdatePolicyBase() {};
  
private:
  //--------------------------- Private Interface ---------------------------//
  FieldUpdatePolicyBase(const FieldUpdatePolicyBase& rhs);
  FieldUpdatePolicyBase& operator=(const FieldUpdatePolicyBase& rhs);
};

}

#ifndef __GCCXML__
#include "FieldUpdatePolicyBaseInline.hh"
#endif

#endif
