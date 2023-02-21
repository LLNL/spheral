text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ManufacturedSolution.cc"

namespace Spheral {
  template class ManufacturedFunction<Dim<%(ndim)s>>;
  template class ManufacturedSteadyStateFunction<Dim<%(ndim)s>>;
  template class ManufacturedConstantFunction<Dim<%(ndim)s>>;
  template class ManufacturedSinusoidalFunction<Dim<%(ndim)s>>;
  template class ManufacturedWaveFunction<Dim<%(ndim)s>>;
  template class ManufacturedTransportSolution<Dim<%(ndim)s>>;
}
"""
