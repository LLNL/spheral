//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Mesh/computeGenerators.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

  template void computeGenerators<Dim<1>,
                                  vector<NodeList<Dim<1> >*>::iterator,
                                  vector<Boundary<Dim<1> >*>::iterator>
  (vector<NodeList<Dim<1> >*>::iterator nodeListBegin,
   vector<NodeList<Dim<1> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
   vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<1>::Vector& xmin,
   const Dim<1>::Vector& xmax,
   vector<Dim<1>::Vector>& positions,
   vector<Dim<1>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<1>,
                                  vector<const NodeList<Dim<1> >*>::iterator,
                                  vector<Boundary<Dim<1> >*>::iterator>
  (vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
   vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
   vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<1>::Vector& xmin,
   const Dim<1>::Vector& xmax,
   vector<Dim<1>::Vector>& positions,
   vector<Dim<1>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<1>,
                                  vector<const NodeList<Dim<1> >*>::iterator,
                                  vector<Boundary<Dim<1> >*>::const_iterator>
  (vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
   vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<1> >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim<1> >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<1>::Vector& xmin,
   const Dim<1>::Vector& xmax,
   vector<Dim<1>::Vector>& positions,
   vector<Dim<1>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<1>,
                                  vector<NodeList<Dim<1> >*>::const_iterator,
                                  vector<Boundary<Dim<1> >*>::const_iterator>
  (vector<NodeList<Dim<1> >*>::const_iterator nodeListBegin,
   vector<NodeList<Dim<1> >*>::const_iterator nodeListEnd,
   vector<Boundary<Dim<1> >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim<1> >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<1>::Vector& xmin,
   const Dim<1>::Vector& xmax,
   vector<Dim<1>::Vector>& positions,
   vector<Dim<1>::SymTensor>& Hs,
   vector<unsigned>& offsets);

#endif

#if defined(SPHERAL_ENABLE_2D)

  template void computeGenerators<Dim<2>,
                                  vector<NodeList<Dim<2> >*>::iterator,
                                  vector<Boundary<Dim<2> >*>::iterator>
  (vector<NodeList<Dim<2> >*>::iterator nodeListBegin,
   vector<NodeList<Dim<2> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
   vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<2>::Vector& xmin,
   const Dim<2>::Vector& xmax,
   vector<Dim<2>::Vector>& positions,
   vector<Dim<2>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<2>,
                                  vector<const NodeList<Dim<2> >*>::iterator,
                                  vector<Boundary<Dim<2> >*>::iterator>
  (vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
   vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
   vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<2>::Vector& xmin,
   const Dim<2>::Vector& xmax,
   vector<Dim<2>::Vector>& positions,
   vector<Dim<2>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<2>,
                                  vector<const NodeList<Dim<2> >*>::iterator,
                                  vector<Boundary<Dim<2> >*>::const_iterator>
  (vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
   vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<2> >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim<2> >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<2>::Vector& xmin,
   const Dim<2>::Vector& xmax,
   vector<Dim<2>::Vector>& positions,
   vector<Dim<2>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<2>,
                                  vector<NodeList<Dim<2> >*>::const_iterator,
                                  vector<Boundary<Dim<2> >*>::const_iterator>
  (vector<NodeList<Dim<2> >*>::const_iterator nodeListBegin,
   vector<NodeList<Dim<2> >*>::const_iterator nodeListEnd,
   vector<Boundary<Dim<2> >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim<2> >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<2>::Vector& xmin,
   const Dim<2>::Vector& xmax,
   vector<Dim<2>::Vector>& positions,
   vector<Dim<2>::SymTensor>& Hs,
   vector<unsigned>& offsets);

#endif

#if defined(SPHERAL_ENABLE_3D)

  template void computeGenerators<Dim<3>,
                                  vector<NodeList<Dim<3> >*>::iterator,
                                  vector<Boundary<Dim<3> >*>::iterator>
  (vector<NodeList<Dim<3> >*>::iterator nodeListBegin,
   vector<NodeList<Dim<3> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
   vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<3>::Vector& xmin,
   const Dim<3>::Vector& xmax,
   vector<Dim<3>::Vector>& positions,
   vector<Dim<3>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<3>,
                                  vector<const NodeList<Dim<3> >*>::iterator,
                                  vector<Boundary<Dim<3> >*>::iterator>
  (vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
   vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
   vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<3>::Vector& xmin,
   const Dim<3>::Vector& xmax,
   vector<Dim<3>::Vector>& positions,
   vector<Dim<3>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<3>,
                                  vector<const NodeList<Dim<3> >*>::iterator,
                                  vector<Boundary<Dim<3> >*>::const_iterator>
  (vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
   vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
   vector<Boundary<Dim<3> >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim<3> >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<3>::Vector& xmin,
   const Dim<3>::Vector& xmax,
   vector<Dim<3>::Vector>& positions,
   vector<Dim<3>::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim<3>,
                                  vector<NodeList<Dim<3> >*>::const_iterator,
                                  vector<Boundary<Dim<3> >*>::const_iterator>
  (vector<NodeList<Dim<3> >*>::const_iterator nodeListBegin,
   vector<NodeList<Dim<3> >*>::const_iterator nodeListEnd,
   vector<Boundary<Dim<3> >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim<3> >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim<3>::Vector& xmin,
   const Dim<3>::Vector& xmax,
   vector<Dim<3>::Vector>& positions,
   vector<Dim<3>::SymTensor>& Hs,
   vector<unsigned>& offsets);

#endif
}