//------------------------------------------------------------------------------
// Templated helper method to construct a python list of objects based
// on a given iterator pair.
//------------------------------------------------------------------------------
#ifndef __Spheral_iteratorsAsList_hh__
#define __Spheral_iteratorsAsList_hh__

#ifndef __GCCXML__
#include "boost/python.hpp"
#else
#include "fakeboost.hh"
#endif

namespace Spheral {

template<typename IteratorType>
inline
boost::python::list
iteratorsAsListByRef(IteratorType begin, IteratorType end) {
  boost::python::list result;
  for (IteratorType itr = begin; itr != end; ++itr) {
    result.append(boost::ref(*itr));
  }
  return result;
}

template<typename IteratorType>
inline
boost::python::list
iteratorsAsListByValue(IteratorType begin, IteratorType end) {
  boost::python::list result;
  for (IteratorType itr = begin; itr != end; ++itr) {
    result.append(*itr);
  }
  return result;
}

}

#endif
