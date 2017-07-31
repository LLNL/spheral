//---------------------------------Spheral++----------------------------------//
// ConstantYVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantYVelocityBoundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"

#include "Utilities/DBC.hh"

namespace Spheral {
namespace BoundarySpace {

using namespace std;

using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using FileIOSpace::FileIO;

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantYVelocityBoundary<Dimension>::
ConstantYVelocityBoundary(const NodeList<Dimension>& nodeList,
                          const vector<int>& nodeIndices):
  ConstantVelocityBoundary<Dimension>(nodeList, nodeIndices) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantYVelocityBoundary<Dimension>::~ConstantYVelocityBoundary() {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantYVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(this->valid());

  // Is this field the velocity on the NodeList we're watching?
  if (field.nodeListPtr() == &(this->nodeList()) &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    int i = 0;
    const vector<int> nodeIDs = this->nodeIndices();
    for (vector<int>::const_iterator itr = nodeIDs.begin();
         itr < nodeIDs.end();
         ++itr, ++i) {
      CHECK(*itr < field.numElements());
      CHECK(i < this->velocityCondition().size());
      field[*itr].y(this->velocityCondition()[i].y());
    }
  }
}

}
}

