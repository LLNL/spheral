text = """
// Define a CPP macro for specializations in the .cc file.
#define SPHERAL%(ndim)sDINSTANTIATION

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/ASPHSmoothingScale.cc"
#include "NodeList/ASPHSmoothingScalev2.cc"

namespace Spheral {
  template class ASPHSmoothingScale<Dim<%(ndim)s>>;
  template class ASPHSmoothingScalev2<Dim<%(ndim)s>>;
}
"""
