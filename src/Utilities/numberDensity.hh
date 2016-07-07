//---------------------------------Spheral++----------------------------------//
// numberDensity
//
// Compute the number density using the ASPH sum of the nodes.
//
// Created by JMO, Tue Feb  9 13:51:33 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_numberDensity__
#define __Spheral_numberDensity__

#include <vector>

// Forward declarations.
namespace Spheral {
namespace FieldSpace {
  template<typename Dimension, typename Value> class FieldList;
}
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}
namespace KernelSpace {
  template<typename Dimension> class TableKernel;
}

//------------------------------------------------------------------------------
// Return the number density for all nodes in the DataBase.
// Note this method assumes that all boundary conditions, ghost nodes, and 
// neighbor information is up to date.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldSpace::FieldList<Dimension, typename Dimension::Scalar>
numberDensity(const DataBaseSpace::DataBase<Dimension>& dataBase,
              const KernelSpace::TableKernel<Dimension>& W);

}

#endif

