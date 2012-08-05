#include "NodeList/SphNodeList.hh"
#include "Boundary/ReflectingBoundary.hh"
#include "Material/CGSUnits.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/GammaLawGas.hh"

#include "Geometry/GeomPlane.cc"
#include "Geometry/GeomVector.cc"
#include "Geometry/GeomTensor.cc"
#include "NodeList/NodeList.cc"
#include "NodeList/FluidNodeList.cc"
#include "NodeList/SphNodeList.cc"
#include "Material/CGSUnits.cc"
#include "Material/MKSUnits.cc"
#include "Material/CosmologicalUnits.cc"
#include "Material/PhysicalConstants.cc"
#include "Material/GammaLawGas.cc"
#include "Neighbor/Neighbor.cc"
#include "Neighbor/NestedGridNeighbor.cc"
#include "Neighbor/GridCellIndex.cc"
#include "Neighbor/GridCellPlane.cc"
#include "Field/FieldBase.cc"
#include "Field/Field.cc"
#include "Boundary/Boundary.cc"
#include "Boundary/PlanarBoundary.cc"
#include "Boundary/ReflectingBoundary.cc"
#include "Kernel/Kernel.cc"
#include "Kernel/BSplineKernel.cc"
#include "Kernel/GaussianKernel.cc"
#include "Kernel/PiGaussianKernel.cc"
#include "Kernel/SuperGaussianKernel.cc"
#include "Kernel/W4SplineKernel.cc"

#include "Geometry/Dimension.cc"

int main() {

  int nxNodes = 10;
  int nyNodes = 10;
  int nzNodes = 10;
  int nNodes = nxNodes*nyNodes*nzNodes;

  double h = 1.0;
  double gamma = 5.0/3.0;
  double mu = 1.0;

  int neighborSearchType = 3; // GatherScatter
  int numGridLevels = 10;
  double topGridCellSize = 2.0;
  Dim<3>::Vector origin(0.0, 0.0, 0.0);
  double kernelExtent = 2.0;

  GammaLawGas<Dim<3>, PhysicalConstants<CGSUnits> > eos(gamma, mu);
  SphNodeList< Dim<3> > nodes(nNodes, eos);

  double dx = 1.0;
  double dy = 1.0;
  double dz = 1.0;

  int nxyNodes = nxNodes*nyNodes;
  Field<Dim<3>, Dim<3>::Vector> &positions = nodes.positions();
  Field<Dim<3>, Dim<3>::Scalar> &Hfield = nodes.Hfield();
  for (int i = 0; i < nzNodes; ++i) {
    for (int j = 0; j < nyNodes; ++j) {
      for (int k = 0; k < nxNodes; ++k) {
        int nodeID = i*nxyNodes + j*nxNodes + k;
        positions(nodeID) = Vector3d((k + 0.5)*dx,
                                     (j + 0.5)*dy,
                                     (i + 0.5)*dz);
        Hfield(nodeID) = h;
      }
    }
  }

  NestedGridNeighbor<Dim<3>, Dim<3>::Scalar>
    nestedNeighbor(nodes,
                   neighborSearchType,
                   numGridLevels,
                   topGridCellSize,
                   origin,
                   kernelExtent);
    

  GeomPlane<Dim<3> > xPlane0(Vector3d(0,0,0),
                                   Vector3d(1,0,0));
  GeomPlane<Dim<3> > xPlane1(Vector3d(10,0,0),
                                   Vector3d(-1,0,0));
  GeomPlane<Dim<3> > yPlane0(Vector3d(0,0,0),
                                   Vector3d(0,1,0));
  GeomPlane<Dim<3> > yPlane1(Vector3d(0,10,0),
                                   Vector3d(0,-1,0));
  GeomPlane<Dim<3> > zPlane0(Vector3d(0,0,0),
                                   Vector3d(0,0,1));
  GeomPlane<Dim<3> > zPlane1(Vector3d(0,0,10),
                                   Vector3d(0,0,-1));

  ReflectingBoundary<Dim<3>, SphNodeList<Dim<3> > > xbc0(xPlane0);
  ReflectingBoundary<Dim<3>, SphNodeList<Dim<3> > > xbc1(xPlane1);
  ReflectingBoundary<Dim<3>, SphNodeList<Dim<3> > > ybc0(yPlane0);
  ReflectingBoundary<Dim<3>, SphNodeList<Dim<3> > > ybc1(yPlane1);
  ReflectingBoundary<Dim<3>, SphNodeList<Dim<3> > > zbc0(zPlane0);
  ReflectingBoundary<Dim<3>, SphNodeList<Dim<3> > > zbc1(zPlane1);

  cerr << "nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes: "
       << nodes.numNodes() << " "
       << nodes.numInternalNodes() << " "
       << nodes.numGhostNodes() << endl;

  cerr << "Applying x0 BC." << endl;
  xbc0.setGhostNodes(nodes);
  cerr << "nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes: "
       << nodes.numNodes() << " "
       << nodes.numInternalNodes() << " "
       << nodes.numGhostNodes() << endl;
  
  cerr << "Applying x1 BC." << endl;
  xbc1.setGhostNodes(nodes);
  cerr << "nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes: "
       << nodes.numNodes() << " "
       << nodes.numInternalNodes() << " "
       << nodes.numGhostNodes() << endl;
  
  cerr << "Applying y0 BC." << endl;
  ybc0.setGhostNodes(nodes);
  cerr << "nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes: "
       << nodes.numNodes() << " "
       << nodes.numInternalNodes() << " "
       << nodes.numGhostNodes() << endl;
  
  cerr << "Applying y1 BC." << endl;
  ybc1.setGhostNodes(nodes);
  cerr << "nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes: "
       << nodes.numNodes() << " "
       << nodes.numInternalNodes() << " "
       << nodes.numGhostNodes() << endl;
  
  cerr << "Applying z0 BC." << endl;
  zbc0.setGhostNodes(nodes);
  cerr << "nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes: "
       << nodes.numNodes() << " "
       << nodes.numInternalNodes() << " "
       << nodes.numGhostNodes() << endl;
  
  cerr << "Applying z1 BC." << endl;
  zbc1.setGhostNodes(nodes);
  cerr << "nodes.numNodes, nodes.numInternalNodes, nodes.numGhostNodes: "
       << nodes.numNodes() << " "
       << nodes.numInternalNodes() << " "
       << nodes.numGhostNodes() << endl;
  
  return 1;
}
