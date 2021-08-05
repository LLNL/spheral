text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Mesh/computeGenerators.cc"

namespace Spheral {

  template void computeGenerators<Dim< %(ndim)s >, 
                                  vector<NeighborNodeList<Dim< %(ndim)s > >*>::iterator,
                                  vector<Boundary<Dim< %(ndim)s > >*>::iterator>
  (vector<NeighborNodeList<Dim< %(ndim)s > >*>::iterator nodeListBegin,
   vector<NeighborNodeList<Dim< %(ndim)s > >*>::iterator nodeListEnd,
   vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryBegin,
   vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim< %(ndim)s >::Vector& xmin,
   const Dim< %(ndim)s >::Vector& xmax,
   vector<Dim< %(ndim)s >::Vector>& positions,
   vector<Dim< %(ndim)s >::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim< %(ndim)s >, 
                                  vector<const NeighborNodeList<Dim< %(ndim)s > >*>::iterator,
                                  vector<Boundary<Dim< %(ndim)s > >*>::iterator>
  (vector<const NeighborNodeList<Dim< %(ndim)s > >*>::iterator nodeListBegin,
   vector<const NeighborNodeList<Dim< %(ndim)s > >*>::iterator nodeListEnd,
   vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryBegin,
   vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim< %(ndim)s >::Vector& xmin,
   const Dim< %(ndim)s >::Vector& xmax,
   vector<Dim< %(ndim)s >::Vector>& positions,
   vector<Dim< %(ndim)s >::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim< %(ndim)s >, 
                                  vector<const NeighborNodeList<Dim< %(ndim)s > >*>::iterator,
                                  vector<Boundary<Dim< %(ndim)s > >*>::const_iterator>
  (vector<const NeighborNodeList<Dim< %(ndim)s > >*>::iterator nodeListBegin,
   vector<const NeighborNodeList<Dim< %(ndim)s > >*>::iterator nodeListEnd,
   vector<Boundary<Dim< %(ndim)s > >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim< %(ndim)s > >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim< %(ndim)s >::Vector& xmin,
   const Dim< %(ndim)s >::Vector& xmax,
   vector<Dim< %(ndim)s >::Vector>& positions,
   vector<Dim< %(ndim)s >::SymTensor>& Hs,
   vector<unsigned>& offsets);

  template void computeGenerators<Dim< %(ndim)s >, 
                                  vector<NeighborNodeList<Dim< %(ndim)s > >*>::const_iterator,
                                  vector<Boundary<Dim< %(ndim)s > >*>::const_iterator>
  (vector<NeighborNodeList<Dim< %(ndim)s > >*>::const_iterator nodeListBegin,
   vector<NeighborNodeList<Dim< %(ndim)s > >*>::const_iterator nodeListEnd,
   vector<Boundary<Dim< %(ndim)s > >*>::const_iterator boundaryBegin,
   vector<Boundary<Dim< %(ndim)s > >*>::const_iterator boundaryEnd,
   const bool meshGhostNodes,
   const Dim< %(ndim)s >::Vector& xmin,
   const Dim< %(ndim)s >::Vector& xmax,
   vector<Dim< %(ndim)s >::Vector>& positions,
   vector<Dim< %(ndim)s >::SymTensor>& Hs,
   vector<unsigned>& offsets);

}
"""
