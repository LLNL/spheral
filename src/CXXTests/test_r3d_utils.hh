//------------------------------------------------------------------------------
// test_r3d_utils
//
// A collection of C++ test functions to exercise the r3d_utils methods.
//
// Created by JMO, Wed Jan 11 15:14:32 PST 2017
//------------------------------------------------------------------------------
#ifndef __Spheral_test_r3d_utils__
#define __Spheral_test_r3d_utils__

#include <string>

namespace Spheral {

//------------------------------------------------------------------------------
// Test converting from r2d_poly <-> polygon.
//------------------------------------------------------------------------------
std::string test_polygon_to_r2d_poly();
std::string test_r2d_poly_to_polygon();

//------------------------------------------------------------------------------
// Test converting from r3d_poly <-> polyhedron.
//------------------------------------------------------------------------------
std::string test_polyhedron_to_r3d_poly();
std::string test_r3d_poly_to_polyhedron();

//------------------------------------------------------------------------------
// Test clipping.
//------------------------------------------------------------------------------
std::string test_clip_polygon();
std::string test_orphan_polygon();
std::string test_clip_polyhedron();

}

#endif
