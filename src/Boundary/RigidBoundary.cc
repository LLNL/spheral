//---------------------------------Spheral++----------------------------------//
// RigidBoundary -- Apply a Reflecting boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

#include "RigidBoundary.hh"

using namespace std;

namespace Spheral {
namespace BoundarySpace {

using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::RigidBoundary():
  PlanarBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given plane.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::
RigidBoundary(const GeomPlane<Dimension>& plane):
  PlanarBoundary<Dimension>(plane, plane) {

  // Once the plane has been set, construct the reflection operator.
  mReflectOperator.Identity();
  mReflectOperator -= 2.0*plane.normal().selfdyad();

  // Once we're done the boundary condition should be in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::~RigidBoundary() {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, just perform a copy.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for scalar fields, just perform a copy.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for Vector fields.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  if ((field.name() == HydroFieldNames::position) || 
      (field.name() == HydroFieldNames::velocity))
  {
     for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
        CHECK(ghostItr < this->ghostEnd(nodeList));
        CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
        CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
        field(*ghostItr) = reflectOperator()*field(*controlItr);
     }
  }
  else
  {
     for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
        CHECK(ghostItr < this->ghostEnd(nodeList));
        CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
        CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
        field(*ghostItr) = field(*controlItr);
     }
  }
}

// Specialization for Tensor fields.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
     CHECK(ghostItr < this->ghostEnd(nodeList));
     CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
     CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
     //     field(*ghostItr) = inverseReflectOperator*field(*controlItr)*reflectOperator();
     field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for third rank tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on the set of nodes in violation of the 
// boundary.
//------------------------------------------------------------------------------
// Specialization for int fields.  A no-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
  REQUIRE(valid());
}

// Specialization for scalar fields.  A no-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(valid());
}

// Specialization for vector fields.  Apply the reflection operator to x and v.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  REQUIRE(valid());

  if ((field.name() == HydroFieldNames::position) || 
      (field.name() == HydroFieldNames::velocity))
  {
     const NodeList<Dimension>& nodeList = field.nodeList();
     for (vector<int>::const_iterator itr = this->violationBegin(nodeList);
           itr < this->violationEnd(nodeList); 
           ++itr) {
        CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
        field(*itr) = reflectOperator()*field(*itr);
     }
  }
}

// Specialization for tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  REQUIRE(valid());
}

// Specialization for tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  REQUIRE(valid());
}

// Specialization for third rank tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  REQUIRE(valid());
}

//------------------------------------------------------------------------------
// Dump state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RigidBoundary<Dimension>::
dumpState(FileIOSpace::FileIO& file,
          const std::string& pathName) const {

  // Call the ancestor class.
  PlanarBoundary<Dimension>::dumpState(file, pathName);

  file.write(reflectOperator(), pathName + "/reflectOperator");
}

//------------------------------------------------------------------------------
// Restore state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RigidBoundary<Dimension>::
restoreState(const FileIOSpace::FileIO& file,
             const std::string& pathName) {

  // Call the ancestor class.
  PlanarBoundary<Dimension>::restoreState(file, pathName);

  file.read(mReflectOperator, pathName + "/reflectOperator");
}

//------------------------------------------------------------------------------
// Test if the reflecting boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
RigidBoundary<Dimension>::valid() const {
  return (reflectOperator().Determinant() != 0.0 &&
          PlanarBoundary<Dimension>::valid());
}

}
}
