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
  RedistributionNotification(Object& object,                             // Main object
                             NullArgMemberFunctionType methodToNotify):  // member function to notify *after* redistribution
    RedistributionNotificationHandle(),
    mObjectPtr(&object),
    mPreMethodToNotify(NULL),
    mMethodToNotify(methodToNotify) {}

  RedistributionNotification(Object& object,                             // Main object
                             NullArgMemberFunctionType preMethodToNotify,// member function to notify *before* redistribution
                             NullArgMemberFunctionType methodToNotify):  // member function to notify *after* redistribution
    RedistributionNotificationHandle(),
    mObjectPtr(&object),
    mPreMethodToNotify(preMethodToNotify),
    mMethodToNotify(methodToNotify) {}

  // Destructor.
  virtual ~RedistributionNotification() {}

  // Required by RedistributionNotificationHandle.
  virtual void notifyBeforeRedistribution() const { if (mPreMethodToNotify != NULL) (mObjectPtr->*(mPreMethodToNotify))(); }
  virtual void notifyAfterRedistribution() const { (mObjectPtr->*(mMethodToNotify))(); }

private:
  //--------------------------- Private Interface ---------------------------//
  Object* mObjectPtr;
  NullArgMemberFunctionType mPreMethodToNotify, mMethodToNotify;
};

}

#endif
