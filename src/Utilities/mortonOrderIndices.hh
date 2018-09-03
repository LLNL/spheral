//---------------------------------Spheral++----------------------------------//
// mortonOrderIndices
//
// Compute the Morton ordered hashed indices for the given set of NodeLists.
// 
// Algorithm described in
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Fri Dec 19 14:58:23 PST 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_mortonOrderIndices__
#define __Spheral_mortonOrderIndices__

#include "Utilities/KeyTraits.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
FieldList<Dimension, typename KeyTraits::Key>
mortonOrderIndices(const FieldList<Dimension, typename Dimension::Vector>& positions);

template<typename Dimension>
FieldList<Dimension, typename KeyTraits::Key>
mortonOrderIndices(const DataBase<Dimension>& dataBase);

// Special version allowing the user to pass a mask indicating nodes
// to ignore in the ordering.
template<typename Dimension>
FieldList<Dimension, typename KeyTraits::Key>
mortonOrderIndices(const DataBase<Dimension>& dataBase,
                   const FieldList<Dimension, int>& mask);

}

#endif

