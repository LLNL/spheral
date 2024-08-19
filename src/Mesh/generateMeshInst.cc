//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Mesh/generateMesh.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

  template void generateMesh<Dim<1>,
                             vector<NodeList<Dim<1> >*>::iterator,
                             vector<Boundary<Dim<1> >*>::iterator>
  (const vector<NodeList<Dim<1> >*>::iterator nodeListBegin,
   const vector<NodeList<Dim<1> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
   const vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
   const Dim<1>::Vector& xmin,
   const Dim<1>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<1> >& mesh,
   NodeList<Dim<1> >& voidNodes);

  template void generateMesh<Dim<1>,
                             vector<const NodeList<Dim<1> >*>::iterator,
                             vector<Boundary<Dim<1> >*>::const_iterator>
  (const vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<1> >*>::const_iterator boundaryBegin,
   const vector<Boundary<Dim<1> >*>::const_iterator boundaryEnd,
   const Dim<1>::Vector& xmin,
   const Dim<1>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<1> >& mesh,
   NodeList<Dim<1> >& voidNodes);

  template void generateMesh<Dim<1>,
                             vector<const NodeList<Dim<1> >*>::iterator,
                             vector<Boundary<Dim<1> >*>::iterator>
  (const vector<const NodeList<Dim<1> >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim<1> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<1> >*>::iterator boundaryBegin,
   const vector<Boundary<Dim<1> >*>::iterator boundaryEnd,
   const Dim<1>::Vector& xmin,
   const Dim<1>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<1> >& mesh,
   NodeList<Dim<1> >& voidNodes);

#endif

#if defined(SPHERAL_ENABLE_2D)

  template void generateMesh<Dim<2>,
                             vector<NodeList<Dim<2> >*>::iterator,
                             vector<Boundary<Dim<2> >*>::iterator>
  (const vector<NodeList<Dim<2> >*>::iterator nodeListBegin,
   const vector<NodeList<Dim<2> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
   const vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
   const Dim<2>::Vector& xmin,
   const Dim<2>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<2> >& mesh,
   NodeList<Dim<2> >& voidNodes);

  template void generateMesh<Dim<2>,
                             vector<const NodeList<Dim<2> >*>::iterator,
                             vector<Boundary<Dim<2> >*>::const_iterator>
  (const vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<2> >*>::const_iterator boundaryBegin,
   const vector<Boundary<Dim<2> >*>::const_iterator boundaryEnd,
   const Dim<2>::Vector& xmin,
   const Dim<2>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<2> >& mesh,
   NodeList<Dim<2> >& voidNodes);

  template void generateMesh<Dim<2>,
                             vector<const NodeList<Dim<2> >*>::iterator,
                             vector<Boundary<Dim<2> >*>::iterator>
  (const vector<const NodeList<Dim<2> >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim<2> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<2> >*>::iterator boundaryBegin,
   const vector<Boundary<Dim<2> >*>::iterator boundaryEnd,
   const Dim<2>::Vector& xmin,
   const Dim<2>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<2> >& mesh,
   NodeList<Dim<2> >& voidNodes);

#endif

#if defined(SPHERAL_ENABLE_3D)

  template void generateMesh<Dim<3>,
                             vector<NodeList<Dim<3> >*>::iterator,
                             vector<Boundary<Dim<3> >*>::iterator>
  (const vector<NodeList<Dim<3> >*>::iterator nodeListBegin,
   const vector<NodeList<Dim<3> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
   const vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
   const Dim<3>::Vector& xmin,
   const Dim<3>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<3> >& mesh,
   NodeList<Dim<3> >& voidNodes);

  template void generateMesh<Dim<3>,
                             vector<const NodeList<Dim<3> >*>::iterator,
                             vector<Boundary<Dim<3> >*>::const_iterator>
  (const vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<3> >*>::const_iterator boundaryBegin,
   const vector<Boundary<Dim<3> >*>::const_iterator boundaryEnd,
   const Dim<3>::Vector& xmin,
   const Dim<3>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<3> >& mesh,
   NodeList<Dim<3> >& voidNodes);

  template void generateMesh<Dim<3>,
                             vector<const NodeList<Dim<3> >*>::iterator,
                             vector<Boundary<Dim<3> >*>::iterator>
  (const vector<const NodeList<Dim<3> >*>::iterator nodeListBegin,
   const vector<const NodeList<Dim<3> >*>::iterator nodeListEnd,
   const vector<Boundary<Dim<3> >*>::iterator boundaryBegin,
   const vector<Boundary<Dim<3> >*>::iterator boundaryEnd,
   const Dim<3>::Vector& xmin,
   const Dim<3>::Vector& xmax,
   const bool meshGhostNodes,
   const bool generateVoid,
   const bool generateParallelConnectivity,
   const bool removeBoundaryZones,
   const double voidThreshold,
   Mesh<Dim<3> >& mesh,
   NodeList<Dim<3> >& voidNodes);

#endif
}