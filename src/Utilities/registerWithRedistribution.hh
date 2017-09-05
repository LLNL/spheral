//---------------------------------Spheral++----------------------------------//
// registerWithRedistribution
// Helper methods to facilitate registering objects with the RedistributionRegistrar.
//
// Created by JMO, Thu May 28 15:46:56 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_registerWithRedistribution__
#define __Spheral_registerWithRedistribution__

#include "RedistributionRegistrar.hh"
#include "RedistributionNotification.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// The data type any client classes must hold onto if they want to register for 
// restarting.
//------------------------------------------------------------------------------
typedef std::shared_ptr<RedistributionNotificationHandle> RedistributionRegistrationType;

//------------------------------------------------------------------------------
// Register the given object for notification of Redistribution.  Such classes
// should store the returned value, since its lifetime determines how long the
// registered object will be notified.
//------------------------------------------------------------------------------
template<typename Object>
RedistributionRegistrationType
registerWithRedistribution(Object& object, 
                           void (Object::* notificationMethod)());

#ifndef __GCCXML__
template<typename Object>
inline
RedistributionRegistrationType
registerWithRedistribution(Object& object, 
                           void (Object::* notificationMethod)()) {
  RedistributionRegistrationType handle(new RedistributionNotification<Object>(object, notificationMethod));
  RedistributionRegistrar& registrar = RedistributionRegistrar::instance();
  registrar.registerRedistributionNotificationHandle(handle);
  return handle;
}
#endif

}

#endif
