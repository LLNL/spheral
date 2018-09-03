//---------------------------------Spheral++----------------------------------//
// NodeListHashMap -- A collection of stuff to simplify using the public 
// my_hash_map class with Spheral NodeList's as keys.
//
// Created by JMO, Sun Oct  3 13:39:17 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeListHashMap__
#define __Spheral_NodeListHashMap__

#include <ext/hash_map>
#include "NodeList/NodeList.hh"

using namespace __gnu_cxx;

// Define a hash function for NodeList pointers.
namespace __gnu_cxx {
  template<>
  struct hash<const Spheral::NodeList< Dim<1> >*> {
    std::size_t operator()(const Spheral::NodeList< Dim<1> >* key) const {
      return (std::size_t) key;
    }
  };

  template<>
  struct hash<const Spheral::NodeList< Dim<2> >*> {
    std::size_t operator()(const Spheral::NodeList< Dim<2> >* key) const {
      return (std::size_t) key;
    }
  };

  template<>
  struct hash<const Spheral::NodeList< Dim<3> >*> {
    std::size_t operator()(const Spheral::NodeList< Dim<3> >* key) const {
      return (std::size_t) key;
    }
  };
}

namespace Spheral {
  // Trait class to define the type of hash_map.
  template<typename Dimension, typename ValueType>
  struct NodeListHashMap {
    typedef hash_map<const NodeList<Dimension>*, ValueType> HashMapType;
  };
}

#endif
