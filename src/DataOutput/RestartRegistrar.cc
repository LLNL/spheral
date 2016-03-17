//---------------------------------Spheral++----------------------------------//
// RestartRegistrar
// A singleton object which maintains a list of the restartable objects in
// Spheral.
//
// Created by JMO, Wed May 27 13:44:32 PDT 2009
//----------------------------------------------------------------------------//

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include "RestartRegistrar.hh"
#include "Utilities/removeElements.hh"

using namespace std;

using Spheral::FileIOSpace::FileIO;

namespace Spheral {
namespace DataOutput {

//------------------------------------------------------------------------------
// Weak pointers don't have operator==, so we have to provide something.
//------------------------------------------------------------------------------
template<typename T>
struct CompareWeakPtr: public binary_function<boost::weak_ptr<T>, boost::weak_ptr<T>, bool> {
  typedef typename binary_function<boost::weak_ptr<T>, boost::weak_ptr<T>, bool>::first_argument_type first_argument_type;
  typedef typename binary_function<boost::weak_ptr<T>, boost::weak_ptr<T>, bool>::second_argument_type second_argument_type;
  typedef typename binary_function<boost::weak_ptr<T>, boost::weak_ptr<T>, bool>::result_type result_type;
  result_type operator()(const boost::weak_ptr<T> lhs,
                         const boost::weak_ptr<T> rhs) const {
    return lhs.lock() == rhs.lock();
  }
};

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
RestartRegistrar&
RestartRegistrar::
instance() {
  return *RestartRegistrar::instancePtr();
}

//------------------------------------------------------------------------------
// Get the instance (pointer).
//------------------------------------------------------------------------------
RestartRegistrar*
RestartRegistrar::
instancePtr() {
   if (mInstancePtr == 0) mInstancePtr = new RestartRegistrar;
   CHECK(mInstancePtr != 0);
   mInstancePtr->removeExpiredPointers();
   return mInstancePtr;
}

//------------------------------------------------------------------------------
// Register a RestartHandle.
//------------------------------------------------------------------------------
void
RestartRegistrar::
registerRestartHandle(boost::shared_ptr<RestartHandle> restartHandlePtr,
                      const unsigned priority) {
  this->removeExpiredPointers();
  CHECK(mPriorities.size() == mRestartHandles.size());
  if (not haveRestartHandle(restartHandlePtr)) {
    priority_iterator itr = upper_bound(mPriorities.begin(), mPriorities.end(), priority);
    const size_t delta = distance(mPriorities.begin(), itr);
    mRestartHandles.insert(mRestartHandles.begin() + delta, restartHandlePtr);
    mPriorities.insert(itr, priority);
  }
  ENSURE(haveRestartHandle(restartHandlePtr));
  ENSURE(mRestartHandles.size() == mPriorities.size());
}

//------------------------------------------------------------------------------
// Unregister a RestartHandle.
//------------------------------------------------------------------------------
void
RestartRegistrar::
unregisterRestartHandle(boost::shared_ptr<RestartHandle> restartHandlePtr) {
  this->removeExpiredPointers();
  VERIFY(haveRestartHandle(restartHandlePtr));
  iterator itr = find(this->begin(), this->end(), restartHandlePtr);
  CHECK(itr != this->end());
  const size_t delta = distance(this->begin(), itr);
  mRestartHandles.erase(itr);
  mPriorities.erase(mPriorities.begin() + delta);
  ENSURE(not haveRestartHandle(restartHandlePtr));
  ENSURE(mRestartHandles.size() == mPriorities.size());
}

//------------------------------------------------------------------------------
// Check whether the given RestartHandle is registered.
//------------------------------------------------------------------------------
bool
RestartRegistrar::
haveRestartHandle(const boost::shared_ptr<RestartHandle> restartHandlePtr) const {
  const_iterator itr = std::find(this->begin(), this->end(), restartHandlePtr);
  return (itr != this->end());
  // const_iterator itr = std::find_if(this->begin(), this->end(), bind2nd(CompareWeakPtr<RestartHandle>(), restartHandlePtr));
  // return (itr != this->end());
  // const_iterator itr = this->begin();
  // while (itr < this->end() and
  //        itr->lock() != restartHandlePtr) ++itr;
  // return (itr != this->end());
}

//------------------------------------------------------------------------------
// Eliminate any pointers to expired objects.
//------------------------------------------------------------------------------
void
RestartRegistrar::
removeExpiredPointers() {
  vector<size_t> expiredIndicies;
  for (size_t i = 0; i != mRestartHandles.size(); ++i) {
    // If there's only one reference to the object, we're it.  Delete it.
    if (mRestartHandles[i].use_count() == 1) expiredIndicies.push_back(i);
  }
  removeElements(mRestartHandles, expiredIndicies);
  removeElements(mPriorities, expiredIndicies);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    ENSURE(mRestartHandles.size() == mPriorities.size());
    for (const_iterator itr = this->begin();
         itr != this->end();
         ++itr ) ENSURE(itr->use_count() > 1);
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Generate unique labels for each of the objects we currently are restarting.
//------------------------------------------------------------------------------
vector<string>
RestartRegistrar::
uniqueLabels() const {
  const_cast<RestartRegistrar*>(this)->removeExpiredPointers();
  vector<string> result;
  unsigned counter = 0;
  for (const_iterator itr = this->begin();
       itr != this->end();
       ++itr) {
    string l = (**itr).label();
    if (find(result.begin(), result.end(), l) != result.end()) {
      stringstream newlabel;
      newlabel << l << "_" << counter;
      l = newlabel.str();
      ++counter;
    }
    CHECK(find(result.begin(), result.end(), l) == result.end());
    result.push_back(l);
  }
  ENSURE(result.size() == mRestartHandles.size());
  return result;
}

//------------------------------------------------------------------------------
// Print out the labels.
//------------------------------------------------------------------------------
void
RestartRegistrar::
printLabels() const {
  const vector<string> labels = this->uniqueLabels();
  for (std::vector<std::string>::const_iterator itr = labels.begin();
       itr != labels.end();
       ++itr) {
    std::cout << *itr << std::endl;
  }
}

//------------------------------------------------------------------------------
// Dump the state of all objects.
//------------------------------------------------------------------------------
void
RestartRegistrar::
dumpState(FileIO& file) const {
  const vector<string> labels = this->uniqueLabels();
  CHECK(labels.size() == mRestartHandles.size());
  for (size_t i = 0; i != labels.size(); ++i) {
    mRestartHandles[i]->dumpState(file, labels[i]);
  }
}

//------------------------------------------------------------------------------
// Restore the state of all objects.
//------------------------------------------------------------------------------
void
RestartRegistrar::
restoreState(const FileIO& file) const {
  const vector<string> labels = this->uniqueLabels();
  CHECK(labels.size() == mRestartHandles.size());
  for (size_t i = 0; i != labels.size(); ++i) {
    mRestartHandles[i]->restoreState(file, labels[i]);
  }
}

//------------------------------------------------------------------------------
// Default constructor (private).
//------------------------------------------------------------------------------
RestartRegistrar::
RestartRegistrar():
   mRestartHandles(),
   mPriorities() {
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
RestartRegistrar::
~RestartRegistrar() {
}

}
}

//------------------------------------------------------------------------------
// Initialize the static instance pointer.
//-----------------------------------------------------------------------------
Spheral::DataOutput::RestartRegistrar* Spheral::DataOutput::RestartRegistrar::mInstancePtr = 0;
