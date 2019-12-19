text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Mesh/generateMesh.cc"

namespace Spheral {

  template void generateMesh<Dim< %(ndim)s >, 
                             vector<NodeList<Dim< %(ndim)s > >*>::iterator,
                             vector<Boundary<Dim< %(ndim)s > >*>::iterator>
  (const vector<NodeList<Dim< %(ndim)s > >*>::iterator nodeListBegin,
   const vector<NodeList<Dim< %(ndim)s > >*>::iterator nodeListEnd,
   const vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryBegin,
   const vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryEnd,
   const Dim< %(ndim)s >::Vector& xmin,
   const Dim< %(ndim)s >::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim< %(ndim)s > >& mesh,
   NodeList<Dim< %(ndim)s > >& voidNodes);

  template void generateMesh<Dim< %(ndim)s >, 
                             vector<const NodeList<Dim< %(ndim)s > >*>::iterator,
                             vector<Boundary<Dim< %(ndim)s > >*>::const_iterator>
  (const vector<const NodeList<Dim< %(ndim)s > >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim< %(ndim)s > >*>::iterator nodeListEnd,
   const vector<Boundary<Dim< %(ndim)s > >*>::const_iterator boundaryBegin,
   const vector<Boundary<Dim< %(ndim)s > >*>::const_iterator boundaryEnd,
   const Dim< %(ndim)s >::Vector& xmin,
   const Dim< %(ndim)s >::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim< %(ndim)s > >& mesh,
   NodeList<Dim< %(ndim)s > >& voidNodes);

  template void generateMesh<Dim< %(ndim)s >, 
                             vector<const NodeList<Dim< %(ndim)s > >*>::iterator,
                             vector<Boundary<Dim< %(ndim)s > >*>::iterator>
  (const vector<const NodeList<Dim< %(ndim)s > >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim< %(ndim)s > >*>::iterator nodeListEnd,
   const vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryBegin,
   const vector<Boundary<Dim< %(ndim)s > >*>::iterator boundaryEnd,
   const Dim< %(ndim)s >::Vector& xmin,
   const Dim< %(ndim)s >::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim< %(ndim)s > >& mesh,
   NodeList<Dim< %(ndim)s > >& voidNodes);

}
"""
