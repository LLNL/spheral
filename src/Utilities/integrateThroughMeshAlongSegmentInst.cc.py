text = """
// Define a CPP macro for specializations in the .cc file.
#define SPHERAL%(ndim)sDINSTANTIATION

//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------
#include "integrateThroughMeshAlongSegment.cc"

namespace Spheral {

  template Dim< %(ndim)s >::Scalar integrateThroughMeshAlongSegment<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const vector<vector<Dim< %(ndim)s >::Scalar> >& values, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const vector<unsigned>& ncells, const Dim< %(ndim)s >::Vector& s0, const Dim< %(ndim)s >::Vector& s1);

}
"""
