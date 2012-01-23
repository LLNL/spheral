#include "ExtendTableKernel.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Wrap the nperhValues parameter as a python list.
//------------------------------------------------------------------------------
boost::python::list
nperhValues1d(const KernelSpace::TableKernel<Spheral::Dim<1> >* self) {
  return Spheral::iteratorsAsListByValue(self->nperhValues().begin(),
                                         self->nperhValues().end());
}

boost::python::list
nperhValues2d(const KernelSpace::TableKernel<Spheral::Dim<2> >* self) {
  return Spheral::iteratorsAsListByValue(self->nperhValues().begin(),
                                         self->nperhValues().end());
}

boost::python::list
nperhValues3d(const KernelSpace::TableKernel<Spheral::Dim<3> >* self) {
  return Spheral::iteratorsAsListByValue(self->nperhValues().begin(),
                                         self->nperhValues().end());
}

//------------------------------------------------------------------------------
// Wrap the WsumValues parameter as a python list.
//------------------------------------------------------------------------------
boost::python::list
WsumValues1d(const KernelSpace::TableKernel<Spheral::Dim<1> >* self) {
  return Spheral::iteratorsAsListByValue(self->WsumValues().begin(),
                                         self->WsumValues().end());
}

boost::python::list
WsumValues2d(const KernelSpace::TableKernel<Spheral::Dim<2> >* self) {
  return Spheral::iteratorsAsListByValue(self->WsumValues().begin(),
                                         self->WsumValues().end());
}

boost::python::list
WsumValues3d(const KernelSpace::TableKernel<Spheral::Dim<3> >* self) {
  return Spheral::iteratorsAsListByValue(self->WsumValues().begin(),
                                         self->WsumValues().end());
}

}
