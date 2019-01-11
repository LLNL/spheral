text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/nthNodalMoment.cc"

namespace Spheral {
  template
  FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> 
  nthNodalMoment<Dim< %(ndim)s >, vector<NodeList<Dim< %(ndim)s > >*>::const_iterator, 0U>(const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListBegin,
                                                                         const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListEnd,
                                                                         const TableKernel<Dim< %(ndim)s > >& W,
                                                                         const bool renormalize);
  template
  FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
  nthNodalMoment<Dim< %(ndim)s >, vector<NodeList<Dim< %(ndim)s > >*>::const_iterator, 1U>(const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListBegin,
                                                                         const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListEnd,
                                                                         const TableKernel<Dim< %(ndim)s > >& W,
                                                                         const bool renormalize);
  template
  FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor> 
  nthNodalMoment<Dim< %(ndim)s >, vector<NodeList<Dim< %(ndim)s > >*>::const_iterator, 2U>(const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListBegin,
                                                                         const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListEnd,
                                                                         const TableKernel<Dim< %(ndim)s > >& W,
                                                                         const bool renormalize);

  template
  void
  zerothAndFirstNodalMoments<Dim< %(ndim)s >, vector<NodeList<Dim< %(ndim)s > >*>::const_iterator>(const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListBegin,
                                                                                 const vector<NodeList<Dim< %(ndim)s > >*>::const_iterator nodeListEnd,
                                                                                 const TableKernel<Dim< %(ndim)s > >& W,
                                                                                 const bool useGradientAsKernel,
                                                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& zerothMoment,
                                                                                 FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& firstMoment);

}
"""
