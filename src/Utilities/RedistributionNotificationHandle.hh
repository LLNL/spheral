//---------------------------------Spheral++----------------------------------//
// RedistributionNotificationHandle
// An untyped interface class to RedistributionNotification<Object>.  This
// allows us to call the registered methods without knowing the type of object 
// we're working with.
//
// Created by JMO, Thu Dec 10 14:20:06 PST 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_RedistributionNotificationHandle__
#define __Spheral_RedistributionNotificationHandle__

namespace Spheral {

class RedistributionNotificationHandle {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  RedistributionNotificationHandle() {}

  // Destructor.
  virtual ~RedistributionNotificationHandle() {}

  //******************************************************************************
  // Methods all notification objects must provide.
  //******************************************************************************
  virtual void notifyBeforeRedistribution() const = 0;
  virtual void notifyAfterRedistribution() const = 0;

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
