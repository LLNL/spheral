//---------------------------------Spheral++----------------------------------//
// RedistributionNotification
// A class templated on the object type
//
// Created by JMO, Thu Dec 10 14:20:06 PST 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_RedistributionNotification__
#define __Spheral_RedistributionNotification__

#include "RedistributionNotificationHandle.hh"

namespace Spheral {

template<typename Object>
class RedistributionNotification: public RedistributionNotificationHandle {
public:
  //--------------------------- Public Interface ---------------------------//
  // Typedefs to define the pointer to member function types.
  typedef void (Object::* NullArgMemberFunctionType)();

  // Constructors.
  RedistributionNotification(Object& object,
                             NullArgMemberFunctionType methodToNotify):
    RedistributionNotificationHandle(),
    mObjectPtr(&object),
    mMethodToNotify(methodToNotify) {}

  // Destructor.
  virtual ~RedistributionNotification() {}

  // Required by RedistributionNotificationHandle.
  virtual void notifyOfRedistribution() const { (mObjectPtr->*(mMethodToNotify))(); }

private:
  //--------------------------- Private Interface ---------------------------//
  Object* mObjectPtr;
  NullArgMemberFunctionType mMethodToNotify;
};

}

#else

// Forward declaration.
namespace Spheral{
  template<typename Object> class RedistributionNotification;
}

#endif
