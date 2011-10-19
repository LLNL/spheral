// BPL includes.
#include "boost/python.hpp"

namespace Spheral {

// Declare the externally compiled extensions.
void wrapSTLVectorPrimitives();
void wrapSTLVectorField();
void wrapSTLVectorFieldList();
void wrapSTLVectorBoundary();
void wrapSTLVectorPhysics();
void wrapSTLVectorDistributed();
void wrapSTDString();
void wrapPairs();

BOOST_PYTHON_MODULE(CXXTypes) {

  wrapSTLVectorPrimitives();
  wrapSTLVectorField();
  wrapSTLVectorFieldList();
  wrapSTLVectorBoundary();
  wrapSTLVectorPhysics();
  wrapSTLVectorDistributed();
  wrapSTDString();
  wrapPairs();

}

}
