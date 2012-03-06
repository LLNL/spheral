#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "BPLWraps/iteratorsAsList.hh"

#ifndef __GCCXML__
#include "BSplineKernel.hh"
#include "W4SplineKernel.hh"
#include "GaussianKernel.hh"
#include "SuperGaussianKernel.hh"
#include "PiGaussianKernel.hh"
#include "SincKernel.hh"
#include "NSincPolynomialKernel.hh"
#include "NBSplineKernel.hh"
#include "HatKernel.hh"
#include "QuarticSplineKernel.hh"
#include "QuinticSplineKernel.hh"
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Wrap the nperhValues parameter as a python list.
//------------------------------------------------------------------------------
boost::python::list
nperhValues1d(const KernelSpace::TableKernel<Spheral::Dim<1> >* self);

boost::python::list
nperhValues2d(const KernelSpace::TableKernel<Spheral::Dim<2> >* self);

boost::python::list
nperhValues3d(const KernelSpace::TableKernel<Spheral::Dim<3> >* self);

//------------------------------------------------------------------------------
// Wrap the WsumValues parameter as a python list.
//------------------------------------------------------------------------------
boost::python::list
WsumValues1d(const KernelSpace::TableKernel<Spheral::Dim<1> >* self);

boost::python::list
WsumValues2d(const KernelSpace::TableKernel<Spheral::Dim<2> >* self);

boost::python::list
WsumValues3d(const KernelSpace::TableKernel<Spheral::Dim<3> >* self);

}
