#include <vector>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "NodeListWrappers.hh"

namespace py = pybind11;
using namespace pybind11::literals;

namespace Spheral {
namespace NodeSpace {

//------------------------------------------------------------------------------
// PyNodeList
//------------------------------------------------------------------------------
template<typename Dimension>
class PyNodeList: public NodeList<Dimension> {
public:
  using NodeList<Dimension>::NodeList;  // inherit constructors

  virtual void deleteNodes(const std::vector<int>& nodeIDs) override {
    PYBIND11_OVERLOAD(void,        // Return type
                      NodeList<Dimension>,    // Parent class
                      deleteNodes, // name of method
                      nodeIDs     // arguments
      );
  }

  virtual std::list<std::vector<char>> packNodeFieldValues(const std::vector<int>& nodeIDs) const override {
    PYBIND11_OVERLOAD(std::list<std::vector<char>>,  // Return type
                      NodeList<Dimension>,                      // Parent class
                      packNodeFieldValues,           // name of method
                      nodeIDs                        // arguments
      );
  }

  virtual void appendInternalNodes(const int numNewNodes, const std::list<std::vector<char>>& packedFieldValues) override {
    PYBIND11_OVERLOAD(void,                        // Return type
                      NodeList<Dimension>,                    // Parent class
                      appendInternalNodes,         // name of method
                      numNewNodes,                 // arguments
                      packedFieldValues
      );
  }

  virtual void reorderNodes(const std::vector<int>& newOrdering) override {
    PYBIND11_OVERLOAD(void,                        // Return type
                      NodeList<Dimension>,                    // Parent class
                      reorderNodes,                // name of method
                      newOrdering                 // arguments
      );
  }

  virtual std::string label() const override {
    PYBIND11_OVERLOAD(std::string,          // Return type
                      NodeList<Dimension>,             // Parent class
                      label                // name of method
      );
  }

  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const override {
    PYBIND11_OVERLOAD(void,          // Return type
                      NodeList<Dimension>,      // Parent class
                      dumpState,     // name of method
                      file,          // arguments
                      pathName
      );
  }

  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName) override {
    PYBIND11_OVERLOAD(void,          // Return type
                      NodeList<Dimension>,             // Parent class
                      restoreState,         // name of method
                      file,                 // arguments
                      pathName
      );
  }

};

}
}

namespace {  // anonymous

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  //............................................................................
  // NodeListRegistrar
  py::class_<Spheral::NodeSpace::NodeList<Dimension>, Spheral::NodeSpace::PyNodeList<Dimension>, std::unique_ptr<Spheral::NodeSpace::NodeList<Dimension>, py::nodelete>>(m, ("NodeList" + suffix).c_str(), py::metaclass())
    ;

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralNodeList) {
  py::module m("SpheralNodeList", "Spheral NodeList module.");

  //............................................................................
  // NodeType
  py::enum_<Spheral::NodeSpace::NodeType>(m, "NodeType")
    .value("InternalNode", Spheral::NodeSpace::InternalNode)
    .value("GhostNode", Spheral::NodeSpace::GhostNode)
    .export_values();

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif

  return m.ptr();
}
