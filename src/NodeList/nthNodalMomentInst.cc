//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "nthNodalMoment.cc"

namespace Spheral {
  namespace NodeSpace {
    // 1-D
    template
    FieldList<Dim<1>, Dim<1>::Scalar> 
    nthNodalMoment<Dim<1>, vector<NodeList<Dim<1> >*>::const_iterator, 0U>(const vector<NodeList<Dim<1> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<1> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<1> >& W,
                                                                           const bool renormalize);
    template
    FieldList<Dim<1>, Dim<1>::Vector> 
    nthNodalMoment<Dim<1>, vector<NodeList<Dim<1> >*>::const_iterator, 1U>(const vector<NodeList<Dim<1> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<1> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<1> >& W,
                                                                           const bool renormalize);
    template
    FieldList<Dim<1>, Dim<1>::SymTensor> 
    nthNodalMoment<Dim<1>, vector<NodeList<Dim<1> >*>::const_iterator, 2U>(const vector<NodeList<Dim<1> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<1> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<1> >& W,
                                                                           const bool renormalize);

    template
    void
    zerothAndFirstNodalMoments<Dim<1>, vector<NodeList<Dim<1> >*>::const_iterator>(const vector<NodeList<Dim<1> >*>::const_iterator nodeListBegin,
                                                                                   const vector<NodeList<Dim<1> >*>::const_iterator nodeListEnd,
                                                                                   const TableKernel<Dim<1> >& W,
                                                                                   const bool useGradientAsKernel,
                                                                                   FieldList<Dim<1>, Dim<1>::Scalar>& zerothMoment,
                                                                                   FieldList<Dim<1>, Dim<1>::Vector>& firstMoment);

    // 2-D
    template
    FieldList<Dim<2>, Dim<2>::Scalar> 
    nthNodalMoment<Dim<2>, vector<NodeList<Dim<2> >*>::const_iterator, 0U>(const vector<NodeList<Dim<2> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<2> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<2> >& W,
                                                                           const bool renormalize);
    template
    FieldList<Dim<2>, Dim<2>::Vector> 
    nthNodalMoment<Dim<2>, vector<NodeList<Dim<2> >*>::const_iterator, 1U>(const vector<NodeList<Dim<2> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<2> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<2> >& W,
                                                                           const bool renormalize);
    template
    FieldList<Dim<2>, Dim<2>::SymTensor> 
    nthNodalMoment<Dim<2>, vector<NodeList<Dim<2> >*>::const_iterator, 2U>(const vector<NodeList<Dim<2> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<2> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<2> >& W,
                                                                           const bool renormalize);

    template
    void
    zerothAndFirstNodalMoments<Dim<2>, vector<NodeList<Dim<2> >*>::const_iterator>(const vector<NodeList<Dim<2> >*>::const_iterator nodeListBegin,
                                                                                   const vector<NodeList<Dim<2> >*>::const_iterator nodeListEnd,
                                                                                   const TableKernel<Dim<2> >& W,
                                                                                   const bool useGradientAsKernel,
                                                                                   FieldList<Dim<2>, Dim<2>::Scalar>& zerothMoment,
                                                                                   FieldList<Dim<2>, Dim<2>::Vector>& firstMoment);

    // 3-D
    template
    FieldList<Dim<3>, Dim<3>::Scalar> 
    nthNodalMoment<Dim<3>, vector<NodeList<Dim<3> >*>::const_iterator, 0U>(const vector<NodeList<Dim<3> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<3> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<3> >& W,
                                                                           const bool renormalize);
    template
    FieldList<Dim<3>, Dim<3>::Vector> 
    nthNodalMoment<Dim<3>, vector<NodeList<Dim<3> >*>::const_iterator, 1U>(const vector<NodeList<Dim<3> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<3> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<3> >& W,
                                                                           const bool renormalize);
    template
    FieldList<Dim<3>, Dim<3>::SymTensor> 
    nthNodalMoment<Dim<3>, vector<NodeList<Dim<3> >*>::const_iterator, 2U>(const vector<NodeList<Dim<3> >*>::const_iterator nodeListBegin,
                                                                           const vector<NodeList<Dim<3> >*>::const_iterator nodeListEnd,
                                                                           const TableKernel<Dim<3> >& W,
                                                                           const bool renormalize);

    template
    void
    zerothAndFirstNodalMoments<Dim<3>, vector<NodeList<Dim<3> >*>::const_iterator>(const vector<NodeList<Dim<3> >*>::const_iterator nodeListBegin,
                                                                                   const vector<NodeList<Dim<3> >*>::const_iterator nodeListEnd,
                                                                                   const TableKernel<Dim<3> >& W,
                                                                                   const bool useGradientAsKernel,
                                                                                   FieldList<Dim<3>, Dim<3>::Scalar>& zerothMoment,
                                                                                   FieldList<Dim<3>, Dim<3>::Vector>& firstMoment);

  }
}
