#include <iostream>
#include "boost/python.hpp"

namespace bp = boost::python;

namespace Spheral {

void booga() {
  std::cout << "Loading Spheralmodules." << std::endl;
}

BOOST_PYTHON_MODULE(Spheralmodules) {
  booga();
  bp::object Geometry = bp::import("Geometry");
  bp::object FileIO = bp::import("FileIO");
  bp::object DataOutput = bp::import("DataOutput");
  bp::object NodeList = bp::import("NodeList");
  bp::object NodeIterators = bp::import("NodeIterators");
  bp::object Field = bp::import("Field");
  bp::object FieldOperations = bp::import("FieldOperations");
  bp::object Kernel = bp::import("Kernel");
  bp::object SplineKernel = bp::import("SplineKernel");
  bp::object Material = bp::import("Material");
  bp::object Neighbor = bp::import("Neighbor");
  bp::object DataBase = bp::import("DataBase");
  bp::object Boundary = bp::import("Boundary");
  bp::object ArtificialViscosity = bp::import("ArtificialViscosity");
  bp::object Physics = bp::import("Physics");
  bp::object Hydro = bp::import("Hydro");
  bp::object ExternalForce = bp::import("ExternalForce");
  bp::object Gravity = bp::import("Gravity");
  bp::object SPHGravity = bp::import("SPHGravity");
  bp::object Integrator = bp::import("Integrator");
  bp::object CXXTypes = bp::import("CXXTypes");
  bp::object Utilities = bp::import("Utilities");

  bp::object main = bp::import("__main__");
  bp::dict maindict = bp::extract<bp::dict>(main.attr("__dict__"));
  maindict.setdefault("Geometry", Geometry);
  maindict.setdefault("CXXTypes", CXXTypes);

}

}
