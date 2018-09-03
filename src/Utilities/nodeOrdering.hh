//---------------------------------Spheral++----------------------------------//
// nodeOrdering
//
// Compute the order that a given set of nodes should be stepped through
// given a FieldList of things to sort them by.
// The FieldList returned is the one to N indexing corresponding to sorting the 
// input in increasing order.
// 
// Created by JMO, Fri Dec 19 16:13:39 PST 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_nodeOrdering__
#define __Spheral_nodeOrdering__

namespace Spheral {

template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension, typename DataType>
FieldList<Dimension, int>
nodeOrdering(const FieldList<Dimension, DataType>& criteria);

}

#endif

