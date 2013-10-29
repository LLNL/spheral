//---------------------------------Spheral++----------------------------------//
// peanoHilbertOrderIndicies
//
// Compute the PeanoHilbert ordered hashed indicies for the given set of NodeLists.
// 
// Algorithm described in
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Sat Dec 20 22:36:58 PST 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_peanoHilbertOrderIndicies__
#define __Spheral_peanoHilbertOrderIndicies__

#include "Utilities/KeyTraits.hh"

namespace Spheral {

// Forward declarations.
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}
namespace FieldSpace {
  template<typename Dimension, typename DataType> class FieldList;
}

template<typename Dimension>
FieldSpace::FieldList<Dimension, typename KeyTraits::Key>
peanoHilbertOrderIndicies(const DataBaseSpace::DataBase<Dimension>& dataBase);

}

#endif

