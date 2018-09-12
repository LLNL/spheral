//---------------------------------Spheral++----------------------------------//
// registerWithRestart
// Helper methods to facilitate registering objects with the RestartRegistrar.
//
// Created by JMO, Thu May 28 15:46:56 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_registerWithRestart__
#define __Spheral_registerWithRestart__

#include "RestartRegistrar.hh"
#include "Restart.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// The data type any client classes must hold onto if they want to register for 
// restarting.
//------------------------------------------------------------------------------
typedef std::shared_ptr<RestartHandle> RestartRegistrationType;

//------------------------------------------------------------------------------
// Register the given object for Restart.  Such classes should store the 
// returned value, since its lifetime determines how long the registered object
// will be counted for restart.
//------------------------------------------------------------------------------
template<typename Object>
inline
RestartRegistrationType
registerWithRestart(Object& object,
                    const unsigned priority = 100) {
  RestartRegistrationType handle(new Restart<Object>(object));
  RestartRegistrar& registrar = RestartRegistrar::instance();
  registrar.registerRestartHandle(handle, priority);
  return handle;
}

}

#endif
