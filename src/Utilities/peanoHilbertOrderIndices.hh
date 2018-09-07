//---------------------------------Spheral++----------------------------------//
// peanoHilbertOrderIndices
//
// Compute the PeanoHilbert ordered hashed indices for the given set of NodeLists.
// 
// Algorithm described in
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Sat Dec 20 22:36:58 PST 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_peanoHilbertOrderIndices__
#define __Spheral_peanoHilbertOrderIndices__

#include "Utilities/KeyTraits.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
FieldList<Dimension, typename KeyTraits::Key>
peanoHilbertOrderIndices(const FieldList<Dimension, typename Dimension::Vector>& positions);

template<typename Dimension>
FieldList<Dimension, typename KeyTraits::Key>
peanoHilbertOrderIndices(const DataBase<Dimension>& dataBase);

}

#endif

