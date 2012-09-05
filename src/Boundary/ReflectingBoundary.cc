//---------------------------------Spheral++----------------------------------//
// ReflectingBoundary -- Apply a Reflecting boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//

#include "ReflectingBoundary.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/innerProduct.hh"
#include "Field/Field.hh"
#include "DBC.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/planarReflectingOperator.hh"

using namespace std;

namespace Spheral {
namespace BoundarySpace {

using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using Geometry::innerProduct;

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ReflectingBoundary<Dimension>::ReflectingBoundary():
  PlanarBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given plane.
//------------------------------------------------------------------------------
template<typename Dimension>
ReflectingBoundary<Dimension>::
ReflectingBoundary(const GeomPlane<Dimension>& plane):
  PlanarBoundary<Dimension>(plane, plane) {

  // Once the plane has been set, construct the reflection operator.
  mReflectOperator = planarReflectingOperator(plane);

  // Once we're done the boundary condition should be in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ReflectingBoundary<Dimension>::~ReflectingBoundary() {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, just perform a copy.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
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
ReflectingBoundary<Dimension>::
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
ReflectingBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

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
    field(*ghostItr) = reflectOperator()*field(*controlItr);
  }
}

// Specialization for Vector3d fields.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {

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
    Dim<3>::Tensor R3(mReflectOperator.xx(), mReflectOperator.xy(), mReflectOperator.xz(),
                      mReflectOperator.yx(), mReflectOperator.yy(), mReflectOperator.yz(),
                      mReflectOperator.zx(), mReflectOperator.zy(), mReflectOperator.zz());
    field(*ghostItr) = R3*field(*controlItr);
  }
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {

  REQUIRE(valid());

  // Duh!  The inverse of a reflection operator is itself.
//   // Get the inverse of the reflection operator, which is just the transpose.
//   Tensor inverseReflectOperator = reflectOperator().Transpose();

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
    field(*ghostItr) = reflectOperator()*(field(*controlItr)*reflectOperator());
  }
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {

  REQUIRE(valid());

  // Duh!  The inverse of a reflection operator is itself.
//   // Get the inverse of the reflection operator, which is just the transpose.
//   Tensor inverseReflectOperator = reflectOperator().Transpose();

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 && *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
//     field(*ghostItr) = (inverseReflectOperator*field(*controlItr)*reflectOperator()).Symmetric();
    field(*ghostItr) = (reflectOperator()*(field(*controlItr)*reflectOperator())).Symmetric();
  }
}

// Specialization for ThirdRankTensor fields.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
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
    field(*ghostItr) = innerProduct<Dimension>(reflectOperator(), innerProduct<Dimension>(field(*controlItr), reflectOperator()));
  }
}

// Specialization for vector<scalar> fields, just perform a copy.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {

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
ReflectingBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
  REQUIRE(valid());
}

// Specialization for scalar fields.  A no-op.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(valid());
}

// Specialization for vector fields.  Apply the reflection operator.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  REQUIRE(valid());

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = this->violationBegin(nodeList);
       itr < this->violationEnd(nodeList); 
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    field(*itr) = reflectOperator()*field(*itr);
  }
}

// Specialization for vector3d fields.  Apply the reflection operator.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector3d>& field) const {
  REQUIRE(valid());

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = this->violationBegin(nodeList);
       itr < this->violationEnd(nodeList); 
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    Dim<3>::Tensor R3(mReflectOperator.xx(), mReflectOperator.xy(), mReflectOperator.xz(),
                      mReflectOperator.yx(), mReflectOperator.yy(), mReflectOperator.yz(),
                      mReflectOperator.zx(), mReflectOperator.zy(), mReflectOperator.zz());
    field(*itr) = R3*field(*itr);
  }
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  REQUIRE(valid());

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = this->violationBegin(nodeList);
       itr < this->violationEnd(nodeList); 
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    field(*itr) = reflectOperator()*field(*itr)*reflectOperator();
  }
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  REQUIRE(valid());

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = this->violationBegin(nodeList);
       itr < this->violationEnd(nodeList); 
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    field(*itr) = (reflectOperator()*field(*itr)*reflectOperator()).Symmetric();
  }
}

// Specialization for third rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  REQUIRE(valid());

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (vector<int>::const_iterator itr = this->violationBegin(nodeList);
       itr < this->violationEnd(nodeList); 
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    field(*itr) = innerProduct<Dimension>(reflectOperator(), innerProduct<Dimension>(field(*itr), reflectOperator()));
  }
}

//------------------------------------------------------------------------------
// Dump state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ReflectingBoundary<Dimension>::
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
ReflectingBoundary<Dimension>::
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
ReflectingBoundary<Dimension>::valid() const {
  return (reflectOperator().Determinant() != 0.0 &&
          PlanarBoundary<Dimension>::valid());
}

}
}
