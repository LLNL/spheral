text = """
//---------------------------------Spheral++----------------------------------//
// WendlandC6Kernel -- .
//
// Created by CDR, Nov 19 2014
//----------------------------------------------------------------------------//
#include "Kernel.hh"
#include "WendlandC6Kernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class WendlandC6Kernel< Dim< %(ndim)s >  >;
}
"""
