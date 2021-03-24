//---------------------------------Spheral++----------------------------------//
// DEMNodeList -- An abstract base class for the NodeLists to represent
//                  fluids.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
//#include "SmoothingScaleBase.hh"
//#include "Material/EquationOfState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/DataBase.hh"
//#include "DataBase/IncrementState.hh"
//#include "DataBase/ReplaceState.hh"
//#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
//#include "Neighbor/ConnectivityMap.hh"
//#include "Utilities/safeInv.hh"
#include "DEMNodeList.hh"

using std::vector;
using std::list;
using std::string;
using std::cerr;
using std::endl;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with given EOS object, along with optional numInternal nodes,
// numGhost nodes, and name.
//------------------------------------------------------------------------------
template<typename Dimension>
DEMNodeList<Dimension>::
DEMNodeList(string name,
              const int numInternal,
              const int numGhost,
              const Scalar hmin,
              const Scalar hmax,
              const Scalar hminratio,
              const Scalar nPerh,
              const int maxNumNeighbors):
  NodeList<Dimension>(name, numInternal, numGhost, hmin, hmax, hminratio, nPerh, maxNumNeighbors),
  mAngularVelocity("angularVelocity", *this){
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DEMNodeList<Dimension>::
~DEMNodeList() {
}

//------------------------------------------------------------------------------
// Set the angular velocity
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
angularVelocity(const Field<Dimension, typename Dimension::Vector>& omega) {
  mAngularVelocity = omega;
  mAngularVelocity.name("angularVelocity");
}


//------------------------------------------------------------------------------
// Calculate and return the volume per node.
//------------------------------------------------------------------------------
//template<typename Dimension>
//void
//DEMNodeList<Dimension>::
//volume(Field<Dimension, typename Dimension::Scalar>& field) const {
//  REQUIRE(field.nodeListPtr() == this);
//  const Field<Dimension, Scalar>& nodeMass = this->mass();
//  const Field<Dimension, Scalar>& nodeDensity = this->angularVelocity();
//  for (size_t i = 0; i != this->numNodes(); ++i) field[i] = nodeMass[i]*safeInv(nodeDensity[i]);
//  field.name(HydroFieldNames::volume);
//}

//------------------------------------------------------------------------------
// Calculate and return the linear momentum per node.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// DEMNodeList<Dimension>::
// linearMomentum(Field<Dimension, typename Dimension::Vector>& field) const {
//   REQUIRE(field.nodeListPtr() == this);
//   const Field<Dimension, Scalar>& nodeMass = this->mass();
//   const Field<Dimension, Vector>& nodeVelocity = this->velocity();
//   for (size_t i = 0; i != this->numNodes(); ++i) field[i] = nodeMass[i]*nodeVelocity[i];
//   field.name(HydroFieldNames::linearMomentum);
// }

//------------------------------------------------------------------------------
// Calculate and return the total energy per node.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// DEMNodeList<Dimension>::
// totalEnergy(Field<Dimension, typename Dimension::Scalar>& field) const {
//   REQUIRE(field.nodeListPtr() == this);
//   const Field<Dimension, Scalar>& nodeMass = this->mass();
//   const Field<Dimension, Vector>& nodeVelocity = this->velocity();
//   const Field<Dimension, Scalar>& nodeEps = this->specificThermalEnergy();
//   for (size_t i = 0; i != this->numNodes(); ++i) 
//     field[i] = nodeMass[i]*(0.5*nodeVelocity[i].magnitude2() + nodeEps[i]);
//   field.name(HydroFieldNames::totalEnergy);
// }

//------------------------------------------------------------------------------
// Dump the current state of the NodeList to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Dump the ancestor class.
  NodeList<Dimension>::dumpState(file, pathName);

  // Dump each of the internal fields of the DEMNodeList.
  file.write(mAngularVelocity, pathName + "/angularVelocity");
}

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMNodeList<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // Restore the ancestor class.
  NodeList<Dimension>::restoreState(file, pathName);

  // Restore each of the internal fields of the DEMNodeList.
  file.read(mAngularVelocity, pathName + "/angularVelocity");
}  

}
