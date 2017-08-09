//---------------------------------Spheral++----------------------------------//
// ConstantRVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

#include "ConstantRVelocityBoundary.hh"

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
ConstantRVelocityBoundary<Dimension>::
ConstantRVelocityBoundary(const NodeList<Dimension>& nodeList,
                          const vector<int>& nodeIndicies):
  ConstantVelocityBoundary<Dimension>(nodeList, nodeIndicies),
  mRadialVelocity() {
  mRadialVelocity.reserve(nodeIndicies.size());
  const Field<Dimension, Vector>& positions = nodeList.positions();
  const Field<Dimension, Vector>& velocities = nodeList.velocity();
  for (vector<int>::const_iterator itr = nodeIndicies.begin();
       itr != nodeIndicies.end();
       ++itr) {
    mRadialVelocity.push_back(velocities(*itr).dot(positions(*itr).unitVector()));
  }
  ENSURE(mRadialVelocity.size() == nodeIndicies.size());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantRVelocityBoundary<Dimension>::~ConstantRVelocityBoundary() {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantRVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(this->valid());

  // Is this field the velocity on the NodeList we're watching?
  if (field.nodeListPtr() == &(this->nodeList()) &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    int j = 0;
    const vector<int> nodeIDs = this->nodeIndices();
    const Field<Dimension, Vector>& positions = this->nodeList().positions();
    CHECK(nodeIDs.size() == mRadialVelocity.size());
    for (vector<int>::const_iterator itr = nodeIDs.begin();
         itr < nodeIDs.end();
         ++itr, ++j) {
      CHECK(*itr < field.numElements());
      CHECK(j < mRadialVelocity.size());
      const int i = *itr;
      const Vector runit = positions(i).unitVector();
      const Vector vperp = field[i] - field[i].dot(runit)*runit;
      field[i] = mRadialVelocity[j]*runit + vperp;
    }
    CHECK(j == nodeIDs.size() and j == mRadialVelocity.size());
  }
}

//------------------------------------------------------------------------------
// Dump the state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantRVelocityBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  ConstantVelocityBoundary<Dimension>::dumpState(file, pathName);
  file.write(mRadialVelocity, pathName + "/radialVelocities");
}

//------------------------------------------------------------------------------
// Read the state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantRVelocityBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  ConstantVelocityBoundary<Dimension>::restoreState(file, pathName);
  file.read(mRadialVelocity, pathName + "/radialVelocities");
}

}
}

