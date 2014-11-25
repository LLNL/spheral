//---------------------------------Spheral++----------------------------------//
// PlanarBoundary -- Abstract base class for the boundary conditions defined
// on planes.  Planar boundaries in general are defined in terms of parallel
// pairs of planes, representing entrance and exit conditions.
//
// Created by JMO, Thu Mar  2 21:34:25 PST 2000
//----------------------------------------------------------------------------//

#include "PlanarBoundary.hh"
#include "mapPositionThroughPlanes.hh"
#include "Geometry/GeomPlane.hh"
#include "NodeList/FluidNodeList.hh"
#include "FileIO/FileIO.hh"
#include "Mesh/Mesh.hh"

#include "Utilities/DBC.hh"

namespace Spheral {
namespace BoundarySpace {

using namespace std;
using NodeSpace::FluidNodeList;
using FileIOSpace::FileIO;
using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using NeighborSpace::Neighbor;
using MeshSpace::Mesh;

//------------------------------------------------------------------------------
// Internal worker method to help with clipping a box range.
//------------------------------------------------------------------------------
// 1D
void
clipBoxWithPlane(const GeomPlane<Dim<1> >& plane,
                 Dim<1>::Vector& point) {
  if (plane.compare(point) == 1) point = plane.point();
}

// 2D
// For now we only check cardinally aligned planes.
void
clipBoxWithPlane(const GeomPlane<Dim<2> >& plane,
                 Dim<2>::Vector& point) {
  typedef Dim<2>::Vector Vector;
  if (plane.compare(point) == 1) {
    if (fuzzyEqual(std::abs(plane.normal().x()), 1.0)) {
      point.x(plane.point().x());
    } else if (fuzzyEqual(std::abs(plane.normal().y()), 1.0)) {
      point.y(plane.point().y());
    }
  }
}

// 3D
// For now we only check cardinally aligned planes.
void
clipBoxWithPlane(const GeomPlane<Dim<3> >& plane,
                 Dim<3>::Vector& point) {
  typedef Dim<3>::Vector Vector;
  if (plane.compare(point) == 1) {
    if (fuzzyEqual(std::abs(plane.normal().x()), 1.0)) {
      point.x(plane.point().x());
    } else if (fuzzyEqual(std::abs(plane.normal().y()), 1.0)) {
      point.y(plane.point().y());
    } else if (fuzzyEqual(std::abs(plane.normal().z()), 1.0)) {
      point.z(plane.point().z());
    }
  }
}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PlanarBoundary<Dimension>::PlanarBoundary():
  Boundary<Dimension>(),
  mEnterPlane(),
  mExitPlane(),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Construct with the given entrance and exit planes.
//------------------------------------------------------------------------------
template<typename Dimension>
PlanarBoundary<Dimension>::
PlanarBoundary(const GeomPlane<Dimension>& enterPlane,
               const GeomPlane<Dimension>& exitPlane):
  Boundary<Dimension>(),
  mEnterPlane(enterPlane),
  mExitPlane(exitPlane),
  mRestart(DataOutput::registerWithRestart(*this)) {
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PlanarBoundary<Dimension>::~PlanarBoundary() {
}

//------------------------------------------------------------------------------
// Determine the set of ghost nodes for the boundary condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::setGhostNodes(NodeList<Dimension>& nodeList) {

  // Remember which node list we are setting the ghost nodes for.
  this->addNodeList(nodeList);
  typename Boundary<Dimension>::BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& controlNodes = boundaryNodes.controlNodes;
  controlNodes = vector<int>();

  // Get the Neighbor object associated with the node list.
  Neighbor<Dimension>& neighbor = nodeList.neighbor();

  // Here we switch between using the Neighbor magic or just doing an O(N)
  // search for anyone who overlaps the exit plane.
  if (true) {

    // Begin by identifying the set of master and neighbor nodes, where master
    // nodes see through the enter plane, and neighbors see through the exit plane.
    neighbor.setMasterList(enterPlane(), exitPlane());

    // Set the list of control nodes.
    // std::copy(neighbor.masterBegin(), neighbor.masterEnd(), std::back_inserter(controlNodes));
    std::copy(neighbor.coarseNeighborBegin(), neighbor.coarseNeighborEnd(), std::back_inserter(controlNodes));

  } else {

    const unsigned n = nodeList.numNodes();
    const double kernelExtent = neighbor.kernelExtent();
    const Field<Dimension, Vector>& pos = nodeList.positions();
    const Field<Dimension, SymTensor>& H = nodeList.Hfield();
    for (unsigned i = 0; i != n; ++i) {
      const Vector& ri = pos(i);
      const SymTensor& Hi = H(i);
      // const GeomPlane<Dimension> enterPlanePrime(Hi*(mEnterPlane.point() - ri),
      //                                            (Hi*mEnterPlane.normal()).unitVector());
      const GeomPlane<Dimension> exitPlanePrime(Hi*(mExitPlane.point() - ri),
                                                (Hi*mExitPlane.normal()).unitVector());
      if (exitPlanePrime.minimumDistance(Vector::zero) <= kernelExtent) controlNodes.push_back(i);
      // cerr << " --> " << i << " " << ri << " " << enterPlanePrime.minimumDistance(Vector::zero) << " " << exitPlanePrime.minimumDistance(Vector::zero) << endl;
    }

  }

  // std::sort(controlNodes.begin(), controlNodes.end());
  // controlNodes.erase(std::unique(controlNodes.begin(), controlNodes.end()), controlNodes.end());

  // Set the ghost node indicies to correspond to these control nodes.
  setGhostNodeIndicies(nodeList);

  // Assign the positions and H's to the new ghost nodes.
  updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// Set the ghost nodes for the given NodeList using the given set of control
// nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList, 
              const vector<int>& presetControlNodes) {

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;

  // Add this NodeList, creating space for control & ghost nodes.
  this->addNodeList(nodeList);

  // Get the Neighbor object associated with the node list.
  Neighbor<Dimension>& neighbor = nodeList.neighbor();

  // Set the list of control nodes.
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& controlNodes = boundaryNodes.controlNodes;
  controlNodes = presetControlNodes;

  // Set the ghost node indicies to correspond to these control nodes.
  setGhostNodeIndicies(nodeList);

  // Assign the masses and positions to the new ghost nodes.
  updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  In this case violation is being "behind" the entrance plane,
// where behind is defined in terms of the plane normal.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::setViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  this->addNodeList(nodeList);
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& vNodes = boundaryNodes.violationNodes;
  vNodes.resize(0);

  // Loop over all the internal nodes in the NodeList, and put any that are 
  // below the enter plane in the list of nodes in violation.
  const Field<Dimension, Vector>& positions = nodeList.positions();
  for (int nodeID = 0; nodeID < nodeList.numInternalNodes(); ++nodeID) {
    if (positions(nodeID) < mEnterPlane) vNodes.push_back(nodeID);
  }

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::updateViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the set of violation nodes for this NodeList.
  const vector<int>& vNodes = this->violationNodes(nodeList);

  // Loop over these nodes, and reset their positions to valid values.
  Field<Dimension, Vector>& positions = nodeList.positions();
  for (vector<int>::const_iterator itr = vNodes.begin();
       itr < vNodes.end();
       ++itr) {
    CHECK(positions(*itr) <= enterPlane());
    positions(*itr) = mapPosition(positions(*itr), mEnterPlane, mExitPlane);
    // CHECK2((positions(*itr) >= enterPlane()) and
    //        (positions(*itr) >= exitPlane()),
    //        "Bad position mapping: " << *itr << " " << nodeList.firstGhostNode() << " " << positions(*itr));
  }

  // Set the Hfield.
  Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
  this->enforceBoundary(Hfield);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (vNodes);
}    

//------------------------------------------------------------------------------
// Function to map a position through the enter to the exit plane.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
PlanarBoundary<Dimension>::
mapPosition(const Vector& position,
	    const GeomPlane<Dimension>& enterPlane,
	    const GeomPlane<Dimension>& exitPlane) const {
  REQUIRE(enterPlane.valid() and exitPlane.valid());
  return mapPositionThroughPlanes(position, enterPlane, exitPlane);
//   REQUIRE(enterPlane.parallel(exitPlane));
//   const Vector deltaEnter = (position - enterPlane.point()).dot(enterPlane.normal())*enterPlane.normal();
//   const Vector deltaExit = (position - exitPlane.point()).dot(exitPlane.normal())*exitPlane.normal();
//   double sign = (position - enterPlane.point()).dot(enterPlane.normal());
//   sign /= fabs(sign) + 1.0e-20;
//   return position - deltaExit - sign*deltaEnter.magnitude()*exitPlane.normal();
}

//------------------------------------------------------------------------------
// Dump state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
dumpState(FileIOSpace::FileIO& file,
          const std::string& pathName) const {
  file.write(enterPlane(), pathName + "/enterPlane");
  file.write(exitPlane(), pathName + "/exitPlane");
}

//------------------------------------------------------------------------------
// Restore state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
restoreState(const FileIOSpace::FileIO& file,
             const std::string& pathName) {
  file.read(mEnterPlane, pathName + "/enterPlane");
  file.read(mExitPlane, pathName + "/exitPlane");
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PlanarBoundary<Dimension>::valid() const {
  return (enterPlane().valid() and
          exitPlane().valid() and
          enterPlane().parallel(exitPlane()));
}

//------------------------------------------------------------------------------
// Internal method to set the ghost node indicies once the controls are set.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::setGhostNodeIndicies(NodeList<Dimension>& nodeList) {

  // Get the sets of control and ghost nodes.
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& controlNodes = boundaryNodes.controlNodes;
  vector<int>& ghostNodes = boundaryNodes.ghostNodes;

  // Create the ghosts of the control nodes.
  int currentNumGhostNodes = nodeList.numGhostNodes();
  int firstNewGhostNode = nodeList.numNodes();
  nodeList.numGhostNodes(currentNumGhostNodes + controlNodes.size());
  CHECK(nodeList.numNodes() == firstNewGhostNode + controlNodes.size());

  // Fill in the pointers to the ghost nodes.
  ghostNodes.resize(controlNodes.size());
  for (int i = 0; i < controlNodes.size(); ++i) {
    CHECK(i >= 0 and i < controlNodes.size());
    CHECK(i >= 0 and i < ghostNodes.size());
    ghostNodes[i] = firstNewGhostNode + i;
    CHECK(ghostNodes[i] >= nodeList.numInternalNodes() and
          ghostNodes[i] < nodeList.numNodes());
  }
}

//------------------------------------------------------------------------------
// Internal method to set the minimal field values necessary to make the 
// ghost nodes valid.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::updateGhostNodes(NodeList<Dimension>& nodeList) {
  REQUIRE(valid());

  // Get the control and ghost node indicies.
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& controlNodes = boundaryNodes.controlNodes;
  vector<int>& ghostNodes = boundaryNodes.ghostNodes;
  CHECK(controlNodes.size() == ghostNodes.size());

  // Loop over the control/ghost node pairs, and set the ghost positions.
  Field<Dimension, Vector>& positions = nodeList.positions();
  typename vector<int>::const_iterator controlItr = controlNodes.begin();
  typename vector<int>::const_iterator ghostItr = ghostNodes.begin();
  for (;controlItr < controlNodes.end(); ++controlItr, ++ghostItr) {
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK2(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes(),
           "Ghost node index out of bounds:  " << *ghostItr << " " << nodeList.firstGhostNode() << " " << nodeList.numNodes());

//     if (!(positions(*controlItr) >= mEnterPlane and
//           positions(*controlItr) >= mExitPlane)) {
//       cerr << "Node out of bounds "
//            << (*controlItr) << " " 
//            << positions(*controlItr) << endl;
//     }
    // CHECK(positions(*controlItr) >= mEnterPlane and
    //       positions(*controlItr) >= mExitPlane);

    positions(*ghostItr) = mapPosition(positions(*controlItr),
                                       mExitPlane,
                                       mEnterPlane);
    // CHECK(positions(*ghostItr) <= mEnterPlane);
  }

  // Set the Hfield.
  Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
  this->applyGhostBoundary(Hfield);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Clip an input box range.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
clip(typename Dimension::Vector& xmin, typename Dimension::Vector& xmax) const {
  clipBoxWithPlane(mEnterPlane, xmin);
  clipBoxWithPlane(mEnterPlane, xmax);
  clipBoxWithPlane(mExitPlane, xmin);
  clipBoxWithPlane(mExitPlane, xmax);
}

//------------------------------------------------------------------------------
// Find the set of tessellation facets on a plane.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<unsigned>
PlanarBoundary<Dimension>::
facesOnPlane(const Mesh<Dimension>& mesh,
             const GeomPlane<Dimension>& plane,
             const Scalar tol) const {

  typedef typename Mesh<Dimension>::Face Face;
  typedef typename Mesh<Dimension>::Node Node;

  vector<unsigned> result;

  // Flag all the nodes that are on the plane to the given tolerance.
  const unsigned numNodes = mesh.numNodes();
  vector<unsigned> flagNodes(numNodes, 0);
  for (unsigned inode = 0; inode != numNodes; ++inode) {
    flagNodes[inode] = (plane.minimumDistance(mesh.node(inode).position()) <= tol ? 1 : 0);
  }

  // Look for any faces with all its nodes in the plane.
  const unsigned numFaces = mesh.numFaces();
  for (unsigned iface = 0; iface != numFaces; ++iface) {
    const Face& face = mesh.face(iface);
    const vector<unsigned> nodeIDs = face.nodeIDs();
    unsigned i = 0;
    while (i < nodeIDs.size() and flagNodes[nodeIDs[i]] == 1) ++i;
    if (i == nodeIDs.size()) result.push_back(iface);
  }

  return result;
}

}
}
