text = """
// Define a CPP macro for specializations in the .cc file.
#define SPHERAL%(ndim)sDINSTANTIATION

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/ASPHSmoothingScale.cc"

namespace Spheral {
  template class ASPHSmoothingScale<Dim<%(ndim)s>>;
}
"""
