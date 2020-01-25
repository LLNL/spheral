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

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Copy control->target values in a Field
//------------------------------------------------------------------------------
template<typename FieldType>
void
copyFieldValues(FieldType& field,
                const vector<int>& control,
                const vector<int>& target) {
  REQUIRE(control.size() == target.size());
  auto controlItr = control.begin();
  auto targetItr = target.begin();
  for (; controlItr < control.end(); ++controlItr, ++targetItr) {
    CHECK(targetItr < target.end());
    CHECK(*controlItr >= 0 && *controlItr < field.numElements());
    CHECK(*targetItr >= 0 && *targetItr < field.numElements());
    field(*targetItr) = field(*controlItr);
  }
}
  
}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::RigidBoundary():
  ReflectingBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given plane.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::
RigidBoundary(const GeomPlane<Dimension>& plane):
  ReflectingBoundary<Dimension>(plane) {
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
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for scalar fields, just perform a copy.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for Vector fields.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  REQUIRE(this->valid());
  if (field.name() == HydroFieldNames::velocity) {
    ReflectingBoundary<Dimension>::applyGhostBoundary(field);
  } else {
    const auto& nodeList = field.nodeList();
    copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
  }
}

// Specialization for Tensor fields.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for third rank tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for fourth rank tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for fifth rank tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
}

// Specialization for FacetedVolume.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  REQUIRE(this->valid());
  const auto& nodeList = field.nodeList();
  copyFieldValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
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
  REQUIRE(this->valid());
}

// Specialization for scalar fields.  A no-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  REQUIRE(this->valid());
}

// Specialization for vector fields.  Apply the reflection operator to x and v.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  REQUIRE(this->valid());
  if (field.name() == HydroFieldNames::velocity) ReflectingBoundary<Dimension>::applyGhostBoundary(field);
}

// Specialization for tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  REQUIRE(this->valid());
}

// Specialization for tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  REQUIRE(this->valid());
}

// Specialization for third rank tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  REQUIRE(this->valid());
}

// Specialization for fourth rank tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  REQUIRE(this->valid());
}

// Specialization for fifth rank tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  REQUIRE(this->valid());
}

// Specialization for FacetedVolume fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  REQUIRE(this->valid());
}

}
