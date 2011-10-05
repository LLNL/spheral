//---------------------------------Spheral++----------------------------------//
// ConstantXVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantXVelocityBoundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"

#include "DBC.hh"
#include "cdebug.hh"

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
ConstantXVelocityBoundary<Dimension>::
ConstantXVelocityBoundary(const NodeList<Dimension>& nodeList,
                          const vector<int>& nodeIndicies):
  ConstantVelocityBoundary<Dimension>(nodeList, nodeIndicies) {
  cdebug << "ConstantXVelocityBoundary::ConstantXVelocityBoundary" << this << endl;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantXVelocityBoundary<Dimension>::~ConstantXVelocityBoundary() {
  cdebug << "ConstantXVelocityBoundary::~ConstantXVelocityBoundary() " << this << endl;
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantXVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  cdebug << "ConstantXVelocityBoundary::enforceBoundary(VectorField) " << this << endl;

  REQUIRE(this->valid());

  // Is this field the velocity on the NodeList we're watching?
  if (field.nodeListPtr() == &(this->nodeList()) &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    int i = 0;
    const vector<int> nodeIDs = this->nodeIndicies();
    for (vector<int>::const_iterator itr = nodeIDs.begin();
         itr < nodeIDs.end();
         ++itr, ++i) {
      CHECK(*itr < field.numElements());
      CHECK(i < this->velocityCondition().size());
      field[*itr].x(this->velocityCondition()[i].x());
    }
  }
}

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace BoundarySpace {
template class ConstantXVelocityBoundary< Dim<1> >;
template class ConstantXVelocityBoundary< Dim<2> >;
template class ConstantXVelocityBoundary< Dim<3> >;
}
}
