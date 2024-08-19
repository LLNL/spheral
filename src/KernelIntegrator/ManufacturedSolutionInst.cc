//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "ManufacturedSolution.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class ManufacturedFunction<Dim<1>>;
  template class ManufacturedSteadyStateFunction<Dim<1>>;
  template class ManufacturedConstantFunction<Dim<1>>;
  template class ManufacturedSinusoidalFunction<Dim<1>>;
  template class ManufacturedWaveFunction<Dim<1>>;
  template class ManufacturedTransportSolution<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class ManufacturedFunction<Dim<2>>;
  template class ManufacturedSteadyStateFunction<Dim<2>>;
  template class ManufacturedConstantFunction<Dim<2>>;
  template class ManufacturedSinusoidalFunction<Dim<2>>;
  template class ManufacturedWaveFunction<Dim<2>>;
  template class ManufacturedTransportSolution<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class ManufacturedFunction<Dim<3>>;
  template class ManufacturedSteadyStateFunction<Dim<3>>;
  template class ManufacturedConstantFunction<Dim<3>>;
  template class ManufacturedSinusoidalFunction<Dim<3>>;
  template class ManufacturedWaveFunction<Dim<3>>;
  template class ManufacturedTransportSolution<Dim<3>>;
#endif
}