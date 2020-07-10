//---------------------------------Spheral++----------------------------------//
// SafeIndexMap.
// This is a version of std::map that provides a const operator[] and 
// disables automatic insertion of elements on indexing.
// 
// When a user tries to index a SafeIndexMap with a Key that is not present,
// a SafeIndexMapKeyError exception is thrown.
//
// Created by JMO, Tue Dec 29 23:06:15 PST 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_SafeIndexMap_hh__
#define __Spheral_SafeIndexMap_hh__

#include <map>
#include <sstream>
#include <string>
#include <exception>

namespace Spheral {

//------------------------------------------------------------------------------
// Exception thrown on key error into SafeIndexMap.
//------------------------------------------------------------------------------
class SafeIndexMapKeyError:
  public std::exception {
public:

  // Constructor
  SafeIndexMapKeyError (const char* text):
    mReason(text) {}

  SafeIndexMapKeyError (const std::string& text):
    mReason(text) {}

  // Destructor
  ~SafeIndexMapKeyError() throw() {}

  // Description of the exception.
  const char* what() const throw() {
    return mReason.c_str();
  }

private:

  // Explanation of the failure.
  std::string mReason;
};

//------------------------------------------------------------------------------
// SafeIndexMap
//------------------------------------------------------------------------------
template<typename Key, typename Value, 
         typename Comparator = std::less<Key> >
class SafeIndexMap:  
  public std::map<Key, Value, Comparator> {

public:

  // Constructors.
  SafeIndexMap();

  // Destructor.
  virtual ~SafeIndexMap();

  // Override the subscript operation.
  Value& operator[](const Key& key);
  const Value& operator[](const Key& key) const;

  // Provide the std::map indexing behaviour with a cumbersome name.
  Value& unsafeIndex(const Key& key);

private:

  // Gripe about an invalid key.
  void mComplainAboutKey(const Key& key) const;
};

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Key, typename Value, typename Comparator>
inline
SafeIndexMap<Key, Value, Comparator>::
SafeIndexMap():
  std::map<Key, Value, Comparator>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Key, typename Value, typename Comparator>
inline
SafeIndexMap<Key, Value, Comparator>::
~SafeIndexMap() {
}

//------------------------------------------------------------------------------
// Subscripting operations.  Complain about nonexisting keys.
//------------------------------------------------------------------------------
template<typename Key, typename Value, typename Comparator>
inline
Value&
SafeIndexMap<Key, Value, Comparator>::
operator[](const Key& key) {
  typename std::map<Key, Value, Comparator>::iterator itr = this->find(key);
  if (itr == this->end()) mComplainAboutKey(key);
  return itr->second;
}

template<typename Key, typename Value, typename Comparator>
inline
const Value&
SafeIndexMap<Key, Value, Comparator>::
operator[](const Key& key) const {
  typename std::map<Key, Value, Comparator>::const_iterator itr = this->find(key);
  if (itr == this->end()) mComplainAboutKey(key);
  return itr->second;
}


//------------------------------------------------------------------------------
// The unsafe normal map subscripting operation.
//------------------------------------------------------------------------------
template<typename Key, typename Value, typename Comparator>
inline
Value&
SafeIndexMap<Key, Value, Comparator>::
unsafeIndex(const Key& key) {
  return std::map<Key, Value, Comparator>::operator[](key);
}

//----------------------------------------------------------------------------
// Complain about an invalid key.
//----------------------------------------------------------------------------
template<typename Key, typename Value, typename Comparator>
inline
void
SafeIndexMap<Key, Value, Comparator>::
mComplainAboutKey(const Key& key) const {
  std::stringstream complaint;
  try {
    complaint << "SafeIndexMap ERROR: requested Key value of "
              << key 
              << " is not in SafeIndexMap "
              << this
              << std::endl;
  } catch (...) {
    complaint << "SafeIndexMap ERROR: requested Key value is not in SafeIndexMap "
              << this
              << std::endl;
  }
  throw SafeIndexMapKeyError(complaint.str());
}

}

#endif
