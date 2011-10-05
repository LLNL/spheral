#include <boost/python.hpp>
#include "Utilities/globalNodeIDs.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/DataBase.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"

using namespace boost::python;
using namespace Spheral;
using namespace Spheral::NodeSpace;
using Spheral::DataBaseSpace::DataBase;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;

namespace Spheral {

//------------------------------------------------------------------------------
// Provide python friendly multiple NodeList versions.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// int
// numGlobalNodes(list listOfNodeLists) {
//   std::vector<const NodeList<Dimension>*> vectorOfNodeLists;
//   const int n = len(listOfNodeLists);
//   for (int i = 0; i != n; ++i) {
//     const NodeList<Dimension>* nodeListPtr = extract<const NodeList<Dimension>*>(listOfNodeLists[i]);
//     vectorOfNodeLists.push_back(nodeListPtr);
//   }
//   return numGlobalNodes(vectorOfNodeLists.begin(),
//                         vectorOfNodeLists.end());
// }

//------------------------------------------------------------------------------
// Wrap the methods for computing global node IDs.
//------------------------------------------------------------------------------
void wrapGlobalNodeIDs() {
  int (*numGlobalNodesNL1d)(const NodeList<Dim<1> >&) = &numGlobalNodes<Dim<1> >;
  int (*numGlobalNodesNL2d)(const NodeList<Dim<2> >&) = &numGlobalNodes<Dim<2> >;
  int (*numGlobalNodesNL3d)(const NodeList<Dim<3> >&) = &numGlobalNodes<Dim<3> >;

  int (*numGlobalNodesDB1d)(const DataBase<Dim<1> >&) = &numGlobalNodes<Dim<1> >;
  int (*numGlobalNodesDB2d)(const DataBase<Dim<2> >&) = &numGlobalNodes<Dim<2> >;
  int (*numGlobalNodesDB3d)(const DataBase<Dim<3> >&) = &numGlobalNodes<Dim<3> >;

  Field<Dim<1>, int> (*globalNodeIDsNL1d)(const NodeList<Dim<1> >&) = &globalNodeIDs<Dim<1> >;
  Field<Dim<2>, int> (*globalNodeIDsNL2d)(const NodeList<Dim<2> >&) = &globalNodeIDs<Dim<2> >;
  Field<Dim<3>, int> (*globalNodeIDsNL3d)(const NodeList<Dim<3> >&) = &globalNodeIDs<Dim<3> >;

  FieldList<Dim<1>, int> (*globalNodeIDsDB1d)(const DataBase<Dim<1> >&) = &globalNodeIDs<Dim<1> >;
  FieldList<Dim<2>, int> (*globalNodeIDsDB2d)(const DataBase<Dim<2> >&) = &globalNodeIDs<Dim<2> >;
  FieldList<Dim<3>, int> (*globalNodeIDsDB3d)(const DataBase<Dim<3> >&) = &globalNodeIDs<Dim<3> >;

  def("numGlobalNodes", numGlobalNodesNL1d, "Find the global number of node in the given NodeList.");
  def("numGlobalNodes", numGlobalNodesNL2d, "Find the global number of node in the given NodeList.");
  def("numGlobalNodes", numGlobalNodesNL3d, "Find the global number of node in the given NodeList.");

  def("numGlobalNodes", numGlobalNodesDB1d, "Find the global number of node in the given DataBase.");
  def("numGlobalNodes", numGlobalNodesDB2d, "Find the global number of node in the given DataBase.");
  def("numGlobalNodes", numGlobalNodesDB3d, "Find the global number of node in the given DataBase.");

  def("globalNodeIDs", globalNodeIDsNL1d, "Compute global node indicies for all nodes in the DataBase.");
  def("globalNodeIDs", globalNodeIDsNL2d, "Compute global node indicies for all nodes in the DataBase.");
  def("globalNodeIDs", globalNodeIDsNL3d, "Compute global node indicies for all nodes in the DataBase.");

  def("globalNodeIDs", globalNodeIDsDB1d, "Compute global node indicies for all nodes in the DataBase.");
  def("globalNodeIDs", globalNodeIDsDB2d, "Compute global node indicies for all nodes in the DataBase.");
  def("globalNodeIDs", globalNodeIDsDB3d, "Compute global node indicies for all nodes in the DataBase.");
}

}
