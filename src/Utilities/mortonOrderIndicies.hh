//---------------------------------Spheral++----------------------------------//
// mortonOrderIndicies
//
// Compute the Morton ordered hashed indicies for the given set of NodeLists.
// 
// Algorithm described in
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Fri Dec 19 14:58:23 PST 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_mortonOrderIndicies__
#define __Spheral_mortonOrderIndicies__

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
mortonOrderIndicies(const DataBaseSpace::DataBase<Dimension>& dataBase);

// Special version allowing the user to pass a mask indicating nodes
// to ignore in the ordering.
template<typename Dimension>
FieldSpace::FieldList<Dimension, typename KeyTraits::Key>
mortonOrderIndicies(const DataBaseSpace::DataBase<Dimension>& dataBase,
                    const FieldSpace::FieldList<Dimension, int>& mask);

}

#endif

