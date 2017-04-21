//---------------------------------Spheral++----------------------------------//
// AxisBoundary -- A specialized Reflecting boundary condition for use along
// the axis in RZ calculations.
// Note: this boundary is automatically constructred by the RZ hydro objects,
//       so the user should not explicitly add this boundary.
//
// Created by JMO, Sun Aug 14 10:20:53 PDT 2016
//----------------------------------------------------------------------------//

#include "AxisBoundaryRZ.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

using namespace std;

namespace Spheral {
namespace BoundarySpace {

using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using MeshSpace::Mesh;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
AxisBoundaryRZ::AxisBoundaryRZ(const double etamin):
  ReflectingBoundary<Dimension>(GeomPlane<Dim<2> >(Dim<2>::Vector(0.0, 0.0),
                                                   Dim<2>::Vector(0.0, 1.0))),
  mEtaMin(etamin) {
  VERIFY2(etamin >= 0.0, "Error: AxisBoundaryRZ requires a minimum eta >= 0.0");
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
AxisBoundaryRZ::~AxisBoundaryRZ() {
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  In this case violation is being "behind" the entrance plane,
// where behind is defined in terms of the plane normal.
//------------------------------------------------------------------------------
void
AxisBoundaryRZ::setViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  typedef Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  this->addNodeList(nodeList);
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& vNodes = boundaryNodes.violationNodes;
  vNodes.resize(0);

  // Loop over all the internal nodes in the NodeList, and put any that are
  // less than etamin in eta space from the axis in violation.
  const Vector runit(0.0, 1.0);
  const Field<Dimension, Vector>& positions = nodeList.positions();
  const Field<Dimension, SymTensor>& H = nodeList.Hfield();
  for (int nodeID = 0; nodeID < nodeList.numInternalNodes(); ++nodeID) {
    const double Hri = (H(nodeID)*runit).y();
    const double etari = Hri*positions(nodeID).y();
    if (etari <= mEtaMin) vNodes.push_back(nodeID);
  }

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
void
AxisBoundaryRZ::updateViolationNodes(NodeList<Dim<2> >& nodeList) {

  // The effective plane we're reflecting from in eta space.
  GeomPlane<Dim<2> > plane(Vector(0.0, mEtaMin), Vector(0.0, 1.0));

  // Get the set of violation nodes for this NodeList.
  const vector<int>& vNodes = this->violationNodes(nodeList);

  // Loop over these nodes, and reset their positions to valid values.
  const Vector runit(0.0, 1.0);
  Field<Dimension, Vector>& positions = nodeList.positions();
  Field<Dimension, SymTensor>& H = nodeList.Hfield();
  for (vector<int>::const_iterator itr = vNodes.begin();
       itr < vNodes.end();
       ++itr) {
    Vector& posi = positions(*itr);
    const double Hri = (H(*itr)*runit).y();
    double etari = Hri*posi.y();
    CHECK(etari <= mEtaMin);
    etari = 2.0*mEtaMin - etari;
    CHECK(etari >= mEtaMin);
    posi.y(etari/Hri);
  }

  // Set the Hfield.
  Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
  this->enforceBoundary(Hfield);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (vNodes);
}    

//------------------------------------------------------------------------------
// Get/set etamin
//------------------------------------------------------------------------------
double AxisBoundaryRZ::etamin() const {
  return mEtaMin;
}

void AxisBoundaryRZ::etamin(const double x) {
  VERIFY2(x >= 0.0, "Error: AxisBoundaryRZ requires a minimum eta >= 0.0");
  mEtaMin = x;
}

}
}
