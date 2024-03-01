//---------------------------------Spheral++----------------------------------//
// HDF5Traits -- A trait class to translate Spheral data types to HDF5 
// counterparts.
//
// Created by JMO, Thu Mar  8 16:30:53 PST 2001
//----------------------------------------------------------------------------//
#ifndef HDF5Traits_HH
#define HDF5Traits_HH

#include "H5Cpp.h"

namespace Spheral {

template<typename DataType> 
struct HDF5Traits {
  // This is the H5 type that the template argument DataType translates to.
  static const CompType Type;
};

}

#endif
