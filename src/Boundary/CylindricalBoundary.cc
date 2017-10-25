//---------------------------------Spheral++----------------------------------//
// CylindricalBoundary -- Create a 3-D boundary around a set of nodes in the
// xy plane to emulate cylindrical (RZ) coordinates.
//
// Created by JMO, Tue Mar 29 13:48:03 PST 2005
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/innerProduct.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/DBC.hh"

#include "CylindricalBoundary.hh"

namespace Spheral {
namespace BoundarySpace {

using namespace std;
using std::max;
using std::min;
using std::abs;

using NodeSpace::FluidNodeList;
using FileIOSpace::FileIO;
using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using Geometry::innerProduct;

//------------------------------------------------------------------------------
// Construct against the given DataBase.
//------------------------------------------------------------------------------
CylindricalBoundary::
CylindricalBoundary(const DataBase<Dim<3> >& dataBase):
  Boundary<Dim<3> >(),
  mDeltaPhi(dataBase.newGlobalFieldList(0.0, "Delta angle for generating ghosts")),
  mGhostPositions(dataBase.newGlobalFieldList(Dim<3>::Vector(), "Ghost node positions")),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
CylindricalBoundary::
~CylindricalBoundary() {
}

//------------------------------------------------------------------------------
// Set the ghost nodes for the given NodeList.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
setGhostNodes(NodeList<Dim<3> >& nodeList) {
  REQUIRE(mDeltaPhi.fieldForNodeList(nodeList) < mDeltaPhi.end());

  // Add this NodeList, creating space for control & ghost nodes.
  addNodeList(nodeList);
  BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  vector<int>& controlNodes = boundaryNodes.controlNodes;
  vector<int>& ghostNodes = boundaryNodes.ghostNodes;
  controlNodes = vector<int>();
  ghostNodes = vector<int>();

  // The fields of ghost positions and delta angles for this NodeList.
  Field<Dim<3>, Scalar>& deltaPhi = **(mDeltaPhi.fieldForNodeList(nodeList));

  // The ghost node indices.
  const int currentNumGhostNodes = nodeList.numGhostNodes();
  int ghostNodeIndex = nodeList.numNodes();

  // Get the state for this NodeList.
  Field<Dimension, Vector>& positions = nodeList.positions();
  Field<Dimension, SymTensor>& H = nodeList.Hfield();

  // Get the number of nodes per h, and the kernel extent.
  const double nPerh = nodeList.nodesPerSmoothingScale();
  const double kernelExtent = nodeList.neighbor().kernelExtent();
  CHECK(distinctlyGreaterThan(nPerh, 0.0));
  CHECK(distinctlyGreaterThan(kernelExtent, 0.0));

  // The vector of ghost positions we'll build up, to assign all at once at the
  // end.
  vector<Vector> ghostPositions;

  // Iterate over the internal nodes.
  for (int i = 0; i != nodeList.numNodes(); ++i) {

    // The state of the node.
    const double ri = positions(i).y();
    const double zi = positions(i).x();
    const SymTensor& Hi = H(i);
    const Scalar hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    CHECK(ri >= 0.0);

    // Have we computed the set of ghost nodes for this node yet?
    if (deltaPhi(i) == 0.0) {
      deltaPhi(i) = angularSpacing(ri, hzi, nPerh, kernelExtent);
    }
    const double dphi = deltaPhi(i);
    CHECK(dphi > 0.0);

    // Create the ghost nodes positions.
    const int nbins = max(1, int(min(M_PI - dphi, 1.5*kernelExtent*hzi/(ri + 1.0e-50))/dphi + 0.5));
    if (!(nbins >= 1 && nbins*dphi < M_PI)) {
      cerr << nbins << " "
           << dphi <<  " "
           << hzi << " "
           << ri << " "
           << endl;
    }
    CHECK(nbins >= 1 && nbins*dphi < M_PI);
    for (int ibin = 0; ibin != nbins; ++ibin) {
      const double phij = (ibin + 1)*dphi;
      controlNodes.push_back(i);
      ghostNodes.push_back(ghostNodeIndex++);
      ghostPositions.push_back(Vector(zi,
                                      ri*cos(phij),
                                      ri*sin(phij)));
      controlNodes.push_back(i);
      ghostNodes.push_back(ghostNodeIndex++);
      ghostPositions.push_back(Vector(zi,
                                      ri*cos(phij),
                                      -ri*sin(phij)));
    }

  }

  // Create the set of ghost nodes on the NodeList.
  const int numNewGhostNodes = ghostNodes.size();
  CHECK(controlNodes.size() == numNewGhostNodes);
  CHECK(ghostPositions.size() == numNewGhostNodes);
  nodeList.numGhostNodes(currentNumGhostNodes + numNewGhostNodes);
  CHECK(ghostNodeIndex == nodeList.numNodes());
  Field<Dimension, Vector>& storedGhostPositions = **(mGhostPositions.fieldForNodeList(nodeList));
  for (int k = 0; k != numNewGhostNodes; ++k) {
    const int i = ghostNodes[k];
    CHECK(i >= nodeList.firstGhostNode() && i < nodeList.numNodes());
    positions(i) = ghostPositions[k];
    storedGhostPositions(i) = ghostPositions[k];
  }

  // Contract checking...
  BEGIN_CONTRACT_SCOPE
  {
    vector<int>::const_iterator controlItr = controlNodes.begin();
    vector<int>::const_iterator ghostItr = ghostNodes.begin();
    for (; controlItr != controlNodes.end(); ++controlItr, ++ghostItr) {
      CHECK(ghostItr < ghostNodes.end());
      const int i = *controlItr;
      const int j = *ghostItr;
      const double ri = positions(i).y();
      const double zi = positions(i).x();
      CHECK(fuzzyEqual(positions(j).x(), zi));
      CHECK(fuzzyEqual(sqrt(FastMath::square(positions(j).y()) + FastMath::square(positions(j).z())), ri));
    }
  }
  END_CONTRACT_SCOPE

  // We can use the normal ghost boundary enforcement to update the
  // H's of the ghost nodes.
  applyGhostBoundary(H);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Update the state of the current ghost nodes (correcting their positions
// and H tensors.)
//------------------------------------------------------------------------------
void
CylindricalBoundary::
updateGhostNodes(NodeList<Dim<3> >& nodeList) {
  REQUIRE(mDeltaPhi.fieldForNodeList(nodeList) < mDeltaPhi.end());

  // Get the control and ghost nodes for this NodeList.
  const BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  const vector<int>& controlNodes = boundaryNodes.controlNodes;
  const vector<int>& ghostNodes = boundaryNodes.ghostNodes;

  // Get the positions for this NodeList.
  Field<Dimension, Vector>& positions = nodeList.positions();
  Field<Dimension, Vector>& ghostPositions = **(mGhostPositions.fieldForNodeList(nodeList));

  // Update the positions, forcing the ghost nodes to exist in a ring of 
  // appropriate radius based on the control nodes.
  CHECK(controlNodes.size() == ghostNodes.size());
  vector<int>::const_iterator controlItr = controlNodes.begin();
  vector<int>::const_iterator ghostItr = ghostNodes.begin();
  for (; controlItr != controlNodes.end(); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostNodes.end());
    const int i = *controlItr;
    const int j = *ghostItr;
    const double ri = positions(i).y();
    const double zi = positions(i).x();
    const double f = ri/sqrt(ghostPositions(j).y()*ghostPositions(j).y() + 
                             ghostPositions(j).z()*ghostPositions(j).z() + 1.0e-50);
    positions(j) = Vector(zi,
                          f*ghostPositions(j).y(),
                          f*ghostPositions(j).z());
    ghostPositions(j) = positions(j);
  }
  CHECK(controlItr == controlNodes.end());
  CHECK(ghostItr == ghostNodes.end());

  // Contract checking...
  BEGIN_CONTRACT_SCOPE
  {
    vector<int>::const_iterator controlItr = controlNodes.begin();
    vector<int>::const_iterator ghostItr = ghostNodes.begin();
    for (; controlItr != controlNodes.end(); ++controlItr, ++ghostItr) {
      CHECK(ghostItr < ghostNodes.end());
      const int i = *controlItr;
      const int j = *ghostItr;
      const double ri = positions(i).y();
      const double zi = positions(i).x();
      CHECK(fuzzyEqual(positions(j).x(), zi));
      CHECK(fuzzyEqual(sqrt(FastMath::square(positions(j).y()) + FastMath::square(positions(j).z())), ri));
    }
  }
  END_CONTRACT_SCOPE

  // We can use the normal ghost boundary enforcement to update the
  // H's of the ghost nodes.
  Field<Dimension, SymTensor>& H = nodeList.Hfield();
  applyGhostBoundary(H);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, just perform a copy.
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, int>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for scalar fields, just perform a copy.
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Scalar>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for Vector fields.
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = reflectOperator(position(*controlItr),
                                       position(*ghostItr))*field(*controlItr);
  }
}

// Specialization for Tensor fields.
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    field(*ghostItr) = T*(field(*controlItr)*T);
  }
}

// Specialization for symmetric tensors.
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    field(*ghostItr) = T*(field(*controlItr)*T).Symmetric();
  }
}

// Specialization for third rank tensors.
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = ghostBegin(nodeList);
  ThirdRankTensor val;
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    val = ThirdRankTensor::zero;
    const ThirdRankTensor& fc = field(*controlItr);
    for (unsigned i = 0; i != Dimension::nDim; ++i) {
      for (unsigned j = 0; j != Dimension::nDim; ++j) {
        for (unsigned k = 0; k != Dimension::nDim; ++k) {
          for (unsigned q = 0; q != Dimension::nDim; ++q) {
            for (unsigned r = 0; r != Dimension::nDim; ++r) {
              for (unsigned s = 0; s != Dimension::nDim; ++s) {
                val(i,j,k) += T(i,q)*T(j,r)*T(k,s)*fc(q,r,s);
              }
            }
          }
        }
      }
    }
    field(*ghostItr) = val;
  }
}

// Specialization for vector<scalar> fields, just perform a copy.
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, std::vector<Dim<3>::Scalar> >& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  In this case we list all internal nodes as being in violation,
// in order to ensure nodes stay rigorously in the xy plane.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
setViolationNodes(NodeList<Dim<3> >& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  vector<int>& vNodes = boundaryNodes.violationNodes;
  vNodes = vector<int>();
  vNodes.reserve(nodeList.numInternalNodes());
  for (int nodeID = 0; nodeID != nodeList.numInternalNodes(); ++nodeID) 
    vNodes.push_back(nodeID);
  CHECK(vNodes.size() == nodeList.numInternalNodes());

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
updateViolationNodes(NodeList<Dim<3> >& nodeList) {

  // Get the set of violation nodes.
  const vector<int>& vNodes = violationNodes(nodeList);

  // State we care about.
  Field<Dim<3>, Vector>& positions = nodeList.positions();
  Field<Dim<3>, Vector>& velocity = nodeList.velocity();
  Field<Dim<3> , SymTensor>& Hfield = nodeList.Hfield();

  // Make sure the positions are sensible.
  for (vector<int>::const_iterator itr = vNodes.begin();
       itr != vNodes.end();
       ++itr) {

    // Has the node crossed to negative r (y)?
    const double ri = positions(*itr).y();
    if (ri < 0.0) {
      positions(*itr).y(-ri);
      velocity(*itr).y(abs(velocity(*itr).y()));
    }
    CHECK(positions(*itr).y() >= 0.0);

    // Make sure there is no component out of the xy plane.
    positions(*itr).z(0.0);

  }
    
  // The normal Field enforcement policies work H.
  this->enforceBoundary(Hfield);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (vNodes);
}    

//------------------------------------------------------------------------------
// Enforce consistency with the boundary condition on fields of different
// DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, no-op.
void
CylindricalBoundary::
enforceBoundary(Field<Dim<3>, int>& field) const {
}

// Specialization for scalar fields, no-op.
void
CylindricalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Scalar>& field) const {
}

// Specialization for Vector fields, forcing Vector's to lie in the xy plane.
void
CylindricalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = violationBegin(nodeList);
       itr != violationEnd(nodeList);
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    field(*itr).z(0.0);
  }
}

// Specialization for Tensor fields, force all off-diagonal components
// to zero.
void
CylindricalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = violationBegin(nodeList);
       itr != violationEnd(nodeList);
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    const double xx = field(*itr).xx();
    const double yy = field(*itr).yy();
    const double zz = field(*itr).zz();
    field(*itr) = Tensor(xx, 0.0, 0.0,
                         0.0, yy, 0.0,
                         0.0, 0.0, zz);
  }
}

// Specialization for SymTensor fields, force all off-diagonal components
// to zero.
void
CylindricalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = violationBegin(nodeList);
       itr != violationEnd(nodeList);
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    const double xx = field(*itr).xx();
    const double yy = field(*itr).yy();
    const double zz = field(*itr).zz();
    field(*itr) = SymTensor(xx, 0.0, 0.0,
                            0.0, yy, 0.0,
                            0.0, 0.0, zz);
  }
}

// Specialization for ThirdRankTensor fields.  Duh?
void
CylindricalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field) const {
}

//------------------------------------------------------------------------------
// Compute the target azimuthal angular spacing.
//------------------------------------------------------------------------------
double
CylindricalBoundary::
angularSpacing(const double ri,
               const double hzi,
               const double nodesPerh,
               const double kernelExtent) {
  REQUIRE(ri >= 0.0);
  REQUIRE(hzi >= 0.0);
  REQUIRE(nodesPerh > 0.0);
  REQUIRE(kernelExtent > 0.0);

  // Choose a minimum effective radius such based on whether we see more than
  // 1/4*Pi around the circumference of the ring.
  const double rmin = kernelExtent*hzi/(0.125*M_PI);
  CHECK(rmin > 0.0);
  const double reff = ri; // max(rmin, ri);

  // Compute the nominal spacing.
  const double delta = hzi/nodesPerh;
  double dphi = 2.0*atan2(delta, 2.0*reff);
  CHECK(dphi > 0.0);

  // Now find the closest value of dphi that will give us an integer + 1/2
  // number of steps in Pi.
  int twoi = max(2, int(2.0*M_PI/dphi));
  if (twoi % 2 != 0) ++twoi;
  CHECK(twoi % 2 == 0);
  const double result = 2.0*M_PI/(twoi + 1);
  ENSURE(result > 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Dump the state to the given file.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mDeltaPhi, pathName + "/deltaPhi");
  file.write(mGhostPositions, pathName + "/ghostPositions");
}

//------------------------------------------------------------------------------
// Read the state from the given file.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mDeltaPhi, pathName + "/deltaPhi");
  file.read(mGhostPositions, pathName + "/ghostPositions");
}

}
}
