//---------------------------------Spheral++----------------------------------//
// RestartRegistrar
// A singleton object which maintains a list of the restartable objects in
// Spheral.
//
// Created by JMO, Wed May 27 13:44:32 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_RestartRegistrar__
#define __Spheral_RestartRegistrar__

#include "RestartHandle.hh"
#include "Utilities/DBC.hh"

#include <vector>
#include <string>
#include <memory> // std::shared_ptr, std::weak_ptr

// Forward declarations.
namespace Spheral {
  class FileIO;
}

namespace Spheral {

class RestartRegistrar {

public:
  //------------------------===== Public Interface =====-----------------------//
  typedef std::vector<std::weak_ptr<RestartHandle> > RestartHandleContainer;
  typedef RestartHandleContainer::const_iterator const_iterator;
  typedef RestartHandleContainer::iterator iterator;

  // Get the instance.
  static RestartRegistrar& instance();
  static RestartRegistrar* instancePtr();

  // Methods for registering a RestartHandle.
  void registerRestartHandle(std::shared_ptr<RestartHandle> restartHandlePtr,
			     const unsigned priority);
  void unregisterRestartHandle(std::shared_ptr<RestartHandle> restartHandlePtr);
  bool haveRestartHandle(const std::shared_ptr<RestartHandle> restartHandlePtr) const;

  // Eliminate any pointers to expired objects.
  void removeExpiredPointers();

  // Iterators over the registered handles.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // Generate unique labels for each of the restartable things.
  std::vector<std::string> uniqueLabels() const;

  // Print out the current ordering and labels for the restart objects.
  void printLabels() const;

  // Cause all registered objects to write their state to the file.
  void dumpState(FileIO& file) const;

  // Cause all registered objects to restore their state from the file.
  void restoreState(const FileIO& file) const;

private:
  //------------------------===== Private Interface =====----------------------//
  typedef std::vector<unsigned> PriorityContainer;
  typedef PriorityContainer::const_iterator const_priority_iterator;
  typedef PriorityContainer::iterator priority_iterator;

  // The one and only instance.
  static RestartRegistrar* mInstancePtr;

  // The list of RestartHandles and their priorities.
  RestartHandleContainer mRestartHandles;
  PriorityContainer mPriorities;

   // No public constructors, destructor, or assignment.
   RestartRegistrar();
   RestartRegistrar(const RestartRegistrar&);
   RestartRegistrar& operator=(const RestartRegistrar&);
   ~RestartRegistrar();
};

}

#include "RestartRegistrarInline.hh"

#else

// Forward declaration.
namespace Spheral {
  class RestartRegistrar;
}

#endif
