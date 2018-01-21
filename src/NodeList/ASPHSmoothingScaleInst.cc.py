text = """
// Define a CPP macro for specializations in the .cc file.
#define SPHERAL%(ndim)sDINSTANTIATION

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ASPHSmoothingScale.cc"

namespace Spheral {
  namespace NodeSpace {
    template class ASPHSmoothingScale<Dim< %(ndim)s > >;
  }
}
"""
