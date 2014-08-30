//---------------------------------Spheral++----------------------------------//
// integrateThroughMeshAlongSegment
//
// Return the result of integrating a quantity along a line segment.
// The quantity here is assumed to be represented a values in a vector<Value>,
// where the vector<Value> is the value of the quantity in a series of cartesian
// cells whose box is defined by by xmin, xmax, and ncells.
//
// We actually pass in a vector<vector<Value> >, which is a progressively refined
// (by factors of 2 in each dimesion) representation of the data.  The idea is that
// we use the finest level with a non-zero value for the value.
//
// Created by JMO, Wed Feb  3 16:03:46 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_integrateThroughMeshAlongSegment__
#define __Spheral_integrateThroughMeshAlongSegment__

#include <vector>

namespace Spheral {

template<typename Dimension, typename Value>
Value
integrateThroughMeshAlongSegment(const std::vector<std::vector<Value> >& values,
                                 const typename Dimension::Vector& xmin,
                                 const typename Dimension::Vector& xmax,
                                 const std::vector<unsigned>& ncells,
                                 const typename Dimension::Vector& s0,
                                 const typename Dimension::Vector& s1);

}

#endif
