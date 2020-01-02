text = """
//---------------------------------Spheral++----------------------------------//
// WendlandC4Kernel -- .
//
// Created by CDR, Nov 5 2014
//----------------------------------------------------------------------------//
#include "Kernel/Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class WendlandC4Kernel< Dim< %(ndim)s >  >;
}
"""
