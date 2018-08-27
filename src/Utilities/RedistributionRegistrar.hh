//---------------------------------Spheral++----------------------------------//
// RedistributionRegistrar
// A singleton object which maintains a list of methods to notify when a 
// redistribution of the nodes has occurred.
//
// Created by JMO, Thu Dec 10 14:07:01 PST 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_RedistributionRegistrar__
#define __Spheral_RedistributionRegistrar__

#include "RedistributionNotificationHandle.hh"
#include "Utilities/DBC.hh"

#include <vector>
#include <memory>

namespace Spheral {

class RedistributionRegistrar {

public:
  //------------------------===== Public Interface =====-----------------------//
  typedef std::vector<std::weak_ptr<RedistributionNotificationHandle> > RedistributionNotificationHandleContainer;
  typedef RedistributionNotificationHandleContainer::const_iterator const_iterator;
  typedef RedistributionNotificationHandleContainer::iterator iterator;

  // Get the instance.
  static RedistributionRegistrar& instance();
  static RedistributionRegistrar* instancePtr();

  // Methods for registering a RedistributionNotificationHandle.
  void registerRedistributionNotificationHandle(std::shared_ptr<RedistributionNotificationHandle> redistributionHandlePtr);
  void unregisterRedistributionNotificationHandle(std::shared_ptr<RedistributionNotificationHandle> redistributionHandlePtr);
  bool haveRedistributionNotificationHandle(std::weak_ptr<RedistributionNotificationHandle> redistributionHandlePtr) const;

  // Eliminate any pointers to expired objects.
  void removeExpiredPointers();

  // Iterators over the registered handles.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // Send out the notifications.
  void preRedistributionNotifications() const;
  void broadcastRedistributionNotifications() const;

private:
  //------------------------===== Private Interface =====----------------------//
  // The one and only instance.
  static RedistributionRegistrar* mInstancePtr;

  // The list of RedistributionNotificationHandles.
  RedistributionNotificationHandleContainer mRedistributionNotificationHandles;

  // No public constructors, destructor, or assignment.
  RedistributionRegistrar();
  RedistributionRegistrar(const RedistributionRegistrar&);
  RedistributionRegistrar& operator=(const RedistributionRegistrar&);
  ~RedistributionRegistrar();
};

}

#ifndef __GCCXML__
#include "RedistributionRegistrarInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  class RedistributionRegistrar;
}

#endif
