//------------------------------------------------------------------------------
// test_r3d_utils
//
// A collection of C++ test functions to exercise the r3d_utils methods.
//
// Created by JMO, Wed Jan 11 15:14:32 PST 2017
//------------------------------------------------------------------------------
#include "test_r3d_utils.hh"
#include "Utilities/r3d_utils.hh"
#include "Geometry/Dimension.hh"

#include <vector>
#include <random>
#include <string>

namespace Spheral {

using std::vector;
using std::string;
using std::to_string;

namespace {  // anonymous

//------------------------------------------------------------------------------
// Return an H (r2d).
//------------------------------------------------------------------------------
r2d_poly construct_H_r2d() {
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Make ouselves an H shaped polygon.
  vector<r2d_rvec2> verts0 = {{0,0},
                              {1,0},
                              {1,1}, 
                              {2,1},
                              {2,0},
                              {3,0},
                              {3,3}, 
                              {2,3}, 
                              {2,2}, 
                              {1,2}, 
                              {1,3}, 
                              {0,3}};
  const unsigned nv0 = verts0.size();
  r2d_poly poly2d;
  r2d_init_poly(&poly2d, &verts0[0], nv0);
  CHECK(r2d_is_good(&poly2d));
  return poly2d;
}

//------------------------------------------------------------------------------
// Return an H (polygon).
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume construct_H_polygon() {
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;
  vector<Vector> vertices0 = {Vector(0,0),
                                      Vector(1,0),
                                      Vector(1,1), 
                                      Vector(2,1),
                                      Vector(2,0),
                                      Vector(3,0),
                                      Vector(3,3), 
                                      Vector(2,3), 
                                      Vector(2,2), 
                                      Vector(1,2), 
                                      Vector(1,3), 
                                      Vector(0,3)};
  const unsigned nv0 = vertices0.size();
  vector<vector<unsigned> > facets0(nv0, vector<unsigned>(2));
  for (unsigned i = 0; i != nv0; ++i) facets0[i] = {i, i+1};
  facets0[nv0-1][1] = 0;
  const FacetedVolume polygon0(vertices0, facets0);
  CHECK(fuzzyEqual(polygon0.volume(), 7.0, 1.0e-10));
  return polygon0;
}
    
//------------------------------------------------------------------------------
// Return a saw-tooth (polygon).
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume construct_saw_polygon() {
    typedef Dim<2>::Vector Vector;
    typedef Dim<2>::FacetedVolume FacetedVolume;
    vector<Vector> vertices0(11);
    for(unsigned i = 0; i<5;++i)
    {
        vertices0[2*i] = Vector(2.0*i,1.0);
        vertices0[2*i+1] = Vector(2.0*i+1.0,2.0);
    }
    vertices0[9] = Vector(8.0,0.0);
    vertices0[10] = Vector(0.0,0.0);
    reverse(vertices0.begin(),vertices0.end());
    const unsigned nv0 = vertices0.size();
    vector<vector<unsigned> > facets0(nv0,vector<unsigned>(2));
    for (unsigned i = 0; i != nv0; ++i) facets0[i] = {i, i+1};
    facets0[nv0-1][1] = 0;
    const FacetedVolume polygon0(vertices0, facets0);
    return polygon0;
}

//------------------------------------------------------------------------------
// Return a pyramid. (r3d)
//------------------------------------------------------------------------------
r3d_poly construct_pyramid_r3d() {
  const unsigned nverts = 5;
  const unsigned nfaces = 5;
  vector<r3d_rvec3> verts = {{0, 0, 0}, 
                             {1, 0, 0},
                             {1, 1, 0},
                             {0, 1, 0},
                             {0.5, 0.5, 1}};
  r3d_int faces[5][4] = {{0, 3, 2, 1},
                         {0, 1, 4},
                         {1, 2, 4},
                         {2, 3, 4},
                         {3, 0, 4}};
  r3d_int** facesp = new r3d_int*[nfaces];
  for (unsigned j = 0; j != nfaces; ++j) {
    const unsigned n = (j == 0 ? 4 : 3);
    facesp[j] = new r3d_int[n];
    for (unsigned k = 0; k != n; ++k) facesp[j][k] = faces[j][k];
  }
  r3d_int nvertsperface[nfaces] = {  // Array of number of vertices per face.
    4, 3, 3, 3, 3
  };
  r3d_poly pyramid3d;
  r3d_init_poly(&pyramid3d, &verts[0], nverts, facesp, nvertsperface, nfaces);
  CHECK(r3d_is_good(&pyramid3d));
  r3d_real vol0;
  r3d_reduce(&pyramid3d, &vol0, 0);
  CHECK2(fuzzyEqual(vol0, 1.0/3.0, 1.0e-10), "Pyramid volume initialization error: " << vol0);
  return pyramid3d;
}

//------------------------------------------------------------------------------
// Return a pyramid. (polyhedron)
//------------------------------------------------------------------------------
Dim<3>::FacetedVolume construct_pyramid_polyhedron() {
  typedef Dim<3>::Vector Vector;
  const unsigned nverts = 5;
  const unsigned nfaces = 5;
  vector<Vector> verts = {Vector(0, 0, 0), 
                          Vector(1, 0, 0),
                          Vector(1, 1, 0),
                          Vector(0, 1, 0),
                          Vector(0.5, 0.5, 1)};
  vector<vector<unsigned> > faces = {{0, 3, 2, 1},
                                     {0, 1, 4},
                                     {1, 2, 4},
                                     {2, 3, 4},
                                     {3, 0, 4}};
  Dim<3>::FacetedVolume pyramid(verts, faces);
  CHECK2(fuzzyEqual(pyramid.volume(), 1.0/3.0, 1.0e-10), "Pyramid volume initialization error: " << pyramid.volume());
  return pyramid;
}

//------------------------------------------------------------------------------
// Return an icosahedron (r3d)
//------------------------------------------------------------------------------
r3d_poly construct_icosahedron_r3d() {
  const unsigned nverts = 12;
  const unsigned nfaces = 20;
  r3d_int faces[nfaces][3] = {
    // 5 faces around point 0
    {0, 11, 5},
    {0, 5, 1},
    {0, 1, 7},
    {0, 7, 10},
    {0, 10, 11},
    // 5 adjacent faces
    {1, 5, 9},
    {5, 11, 4},
    {11, 10, 2},
    {10, 7, 6},
    {7, 1, 8},
    // 5 faces around point 3
    {3, 9, 4},
    {3, 4, 2},
    {3, 2, 6},
    {3, 6, 8},
    {3, 8, 9},
    // 5 adjacent faces
    {4, 9, 5},
    {2, 4, 11},
    {6, 2, 10},
    {8, 6, 7},
    {9, 8, 1},
  };
  r3d_int** facesp = new r3d_int*[nfaces];
  for (unsigned j = 0; j != nfaces; ++j) {
    facesp[j] = new r3d_int[3];
    for (unsigned k = 0; k != 3; ++k) facesp[j][k] = faces[j][k];
  }
  r3d_int nvertsperface[nfaces] = {  // Array of number of vertices per face.
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
  };
  const double t = (1.0 + sqrt(5.0)) / 2.0;
  vector<r3d_rvec3> verts = {           // Array of vertex coordinates.
    {-1,  t,   0},
    { 1,  t,   0},
    {-1, -t,   0},
    { 1, -t,   0},
    { 0, -1,   t},
    { 0,  1,   t},
    { 0, -1,  -t},
    { 0,  1,  -t},
    { t,  0,  -1},
    { t,  0,   1},
    {-t,  0,  -1},
    {-t,  0,   1}
  };
  r3d_poly ico3d;
  r3d_init_poly(&ico3d, &verts[0], nverts, facesp, nvertsperface, nfaces);
  CHECK(r3d_is_good(&ico3d));
  return ico3d;
}

//------------------------------------------------------------------------------
// Return an icosahedron. (polyhedron)
//------------------------------------------------------------------------------
Dim<3>::FacetedVolume construct_icosahedron_polyhedron() {
  typedef Dim<3>::Vector Vector;
  const unsigned nverts = 12;
  const unsigned nfaces = 20;
  vector<vector<unsigned> > facets = {
    // 5 faces around point 0
    {0, 11, 5},
    {0, 5, 1},
    {0, 1, 7},
    {0, 7, 10},
    {0, 10, 11},
    // 5 adjacent faces
    {1, 5, 9},
    {5, 11, 4},
    {11, 10, 2},
    {10, 7, 6},
    {7, 1, 8},
    // 5 faces around point 3
    {3, 9, 4},
    {3, 4, 2},
    {3, 2, 6},
    {3, 6, 8},
    {3, 8, 9},
    // 5 adjacent faces
    {4, 9, 5},
    {2, 4, 11},
    {6, 2, 10},
    {8, 6, 7},
    {9, 8, 1},
  };
  const double t = (1.0 + sqrt(5.0)) / 2.0;
  vector<Vector> verts = {           // Array of vertex coordinates.
    Vector(-1,  t,   0),
    Vector( 1,  t,   0),
    Vector(-1, -t,   0),
    Vector( 1, -t,   0),
    Vector( 0, -1,   t),
    Vector( 0,  1,   t),
    Vector( 0, -1,  -t),
    Vector( 0,  1,  -t),
    Vector( t,  0,  -1),
    Vector( t,  0,   1),
    Vector(-t,  0,  -1),
    Vector(-t,  0,   1)
  };
  return Dim<3>::FacetedVolume(verts, facets);
}

}            // anonymous

//------------------------------------------------------------------------------
// Test converting from polygon -> r2d_poly.
//------------------------------------------------------------------------------
std::string test_polygon_to_r2d_poly() {

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Make ouselves an H shaped polygon.
  const FacetedVolume polygon0 = construct_H_polygon();

  // Convert to a r2d_poly.
  r2d_poly poly2d;
  polygon_to_r2d_poly(polygon0, poly2d);

  // Is it correct?
  if (not r2d_is_good(&poly2d)) return "ERROR: r2d_is_good fails";
  r2d_real area;
  r2d_reduce(&poly2d, &area, 0);
  if (not fuzzyEqual(area, 7.0, 1.0e-10)) return "ERROR: area mismatch: " + to_string(area) + " != 7.0";

  // Must be OK.
  return "OK";
}
    
//------------------------------------------------------------------------------
// Test converting from r2d_poly -> polygon.
//------------------------------------------------------------------------------
std::string test_r2d_poly_to_polygon() {

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Make ouselves an H shaped polygon.
  r2d_poly poly2d = construct_H_r2d();

  // Convert to a polygon.
  Dim<2>::FacetedVolume polygon;
  r2d_poly_to_polygon(poly2d, 1.0e-8, polygon);

  // Is it correct?
  if (not fuzzyEqual(polygon.volume(), 7.0, 1.0e-10)) return "ERROR: area mismatch: " + to_string(polygon.volume()) + " != 7.0";

  // Must be OK.
  return "OK";
}

//------------------------------------------------------------------------------
// Test converting from polyhedron -> r3d_poly.
//------------------------------------------------------------------------------
std::string test_polyhedron_to_r3d_poly() {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  // Cube test.
  {
    vector<Vector> vertices0 = {Vector(1, 1, 1), 
                                Vector(1, 2, 1),
                                Vector(2, 1, 1),
                                Vector(2, 2, 1),
                                Vector(1, 1, 2), 
                                Vector(1, 2, 2),
                                Vector(2, 1, 2),
                                Vector(2, 2, 2)};
    const FacetedVolume cube0(vertices0);   // Builds the convex hull.
    CHECK(fuzzyEqual(cube0.volume(), 1.0, 1.0e-10));

    // Convert to a r3d_poly.
    r3d_poly cube3d;
    polyhedron_to_r3d_poly(cube0, cube3d);

    // Is it correct?
    if (not r3d_is_good(&cube3d)) {
      return "ERROR: r3d_is_good fails for cube";
    }
    r3d_real vol;
    r3d_reduce(&cube3d, &vol, 0);
    if (not fuzzyEqual(vol, 1.0, 1.0e-10)) return "ERROR: volume mismatch for cube: " + to_string(vol) + " != 1.0";
  }

  // Pyramid test.
  {
    const FacetedVolume pyramid = construct_pyramid_polyhedron();
    r3d_poly pyramid3d;
    polyhedron_to_r3d_poly(pyramid, pyramid3d);

    // Is it correct?
    if (not r3d_is_good(&pyramid3d)) {
      return "ERROR: r3d_is_good fails for cube";
    }
    r3d_real vol;
    r3d_reduce(&pyramid3d, &vol, 0);
    if (not fuzzyEqual(vol, 1.0/3.0, 1.0e-10)) return "ERROR: volume mismatch for pyramid: " + to_string(vol) + " != 1.0/3.0";
  }

  // Icosahedron test.
  {
    const FacetedVolume ico0 = construct_icosahedron_polyhedron();

    // Convert to a r3d_poly.
    r3d_poly ico3d;
    polyhedron_to_r3d_poly(ico0, ico3d);

    // Is it correct?
    if (not r3d_is_good(&ico3d)) {
      return "ERROR: r3d_is_good fails for icosahedron";
    }
    const double vol0 = ico0.volume();
    r3d_real vol;
    r3d_reduce(&ico3d, &vol, 0);
    if (not fuzzyEqual(vol, vol0, 1.0e-10)) return "ERROR: volume mismatch for icosahedron: " + to_string(vol) + " != " + to_string(vol0);
  }

  return "OK";
}

//------------------------------------------------------------------------------
// Test converting from r3d_poly -> polyhedron.
//------------------------------------------------------------------------------
std::string test_r3d_poly_to_polyhedron() {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef Dim<3>::FacetedVolume::Facet Facet;

  // Cube test.
  {
    vector<r3d_rvec3> vertices0 = {{0, 0, 0}, 
                                   {1, 1, 1}};
    r3d_poly cube3d;
    r3d_init_box(&cube3d, &vertices0[0]);
    CHECK(r3d_is_good(&cube3d));
    r3d_real vol0;
    r3d_reduce(&cube3d, &vol0, 0);
    CHECK2(fuzzyEqual(vol0, 1.0, 1.0e-10), "Cube volume initialization error: " << vol0);

    FacetedVolume cube;
    r3d_poly_to_polyhedron(cube3d, 1.0e-8, cube);

    // Is it correct?
    if (not fuzzyEqual(cube.volume(), 1.0, 1.0e-10)) return "ERROR: volume mismatch for cube: " + to_string(cube.volume()) + " != 1.0";
  }

  // Pyramid test.
  {
    r3d_poly pyramid3d = construct_pyramid_r3d();
    FacetedVolume pyramid;
    r3d_poly_to_polyhedron(pyramid3d, 1.0e-8, pyramid);

    // Is it correct?
    if (not r3d_is_good(&pyramid3d)) {
      return "ERROR: r3d_is_good fails for pyramid";
    }
    const double vol0 = pyramid.volume();
    r3d_real vol;
    r3d_reduce(&pyramid3d, &vol, 0);
    if (not fuzzyEqual(vol, vol0, 1.0e-10)) return "ERROR: volume mismatch for pyramid: " + to_string(vol) + " != " + to_string(vol0);
  }

  // Icosahedron test.
  {
    r3d_poly ico3d = construct_icosahedron_r3d();

    // Convert to a polyhedron.
    FacetedVolume ico;
    r3d_poly_to_polyhedron(ico3d, 1.0e-12, ico);

    // Is it correct?
    r3d_real vol0;
    r3d_reduce(&ico3d, &vol0, 0);
    if (not fuzzyEqual(ico.volume(), vol0, 1.0e-10)) return "ERROR: volume mismatch for icosahedron: " + to_string(ico.volume()) + " != " + to_string(vol0);
  }

  return "OK";

}

//------------------------------------------------------------------------------
// Test clipping a polygon.
//------------------------------------------------------------------------------
std::string test_clip_polygon() {

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef GeomPlane<Dim<2> > Plane;

  // Make ouselves an H shaped polygon.
  const FacetedVolume H = construct_H_polygon();

  // Make some clip planes.
  vector<Plane> planes = {Plane(Vector(1.5, 1.5), Vector(1, 1).unitVector()),
                          Plane(Vector(0.0, 1.0), Vector(0, 1))};

  // Clip it!
  const FacetedVolume Hclip = clipFacetedVolume(H, planes);

  // Is it correct?
  if (not fuzzyEqual(Hclip.volume(), 3.0, 1.0e-10)) return "ERROR: clipped area mismatch: " + to_string(Hclip.volume()) + " != 3.0";

  // Must be OK.
  return "OK";
}
    
//------------------------------------------------------------------------------
// Test clipping the saw-tooth and orphaning the blades
//------------------------------------------------------------------------------
std::string test_orphan_polygon() {
    typedef Dim<2>::Vector Vector;
    typedef Dim<2>::FacetedVolume FacetedVolume;
    typedef GeomPlane<Dim<2> > Plane;
    
    // Make a saw shaped polygon.
    const FacetedVolume S = construct_saw_polygon();
    
    // One clip plane to orphan the blades
    vector<Plane> planes = {Plane(Vector(0.0,1.5),Vector(0,1))};
    
    // Now clip
    const FacetedVolume Sclip = clipFacetedVolume(S,planes);
    
    // Check correctness
    if (not fuzzyEqual(Sclip.volume(),1.0,1.0e-10)) return "ERROR: orphaned area mismatch: " + to_string(Sclip.volume()) + " != 1.0";
    
    // return
    return "OK";
}

//------------------------------------------------------------------------------
// Test clipping a polyhedron.
//------------------------------------------------------------------------------
std::string test_clip_polyhedron() {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef GeomPlane<Dim<3> > Plane;

  // The starting icosahedron.
  const FacetedVolume ico = construct_icosahedron_polyhedron();

  // Random number generator.
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  // Build our random clip planes.
  vector<Plane> clipPlanes;
  const Vector lengths = ico.xmax() - ico.xmin();
  for (unsigned i = 0; i != 3; ++i) {
    const Vector p = Vector(ico.xmin() + Vector(lengths(0) * distribution(generator),
                                                lengths(1) * distribution(generator),
                                                lengths(2) * distribution(generator)));
    const Vector nhat = -p.unitVector();
    clipPlanes.push_back(Plane(p, nhat));
  }

  // First figure out the expected volume by clipping an r3d icosahedron.
  r3d_real vol0;
  {
    r3d_poly ico3d = construct_icosahedron_r3d();
    r3d_plane planes[3];
    for (unsigned i = 0; i != 3; ++i) {
      planes[i].n.x = clipPlanes[i].normal().x();
      planes[i].n.y = clipPlanes[i].normal().y();
      planes[i].n.z = clipPlanes[i].normal().z();
      planes[i].d = -clipPlanes[i].normal().dot(clipPlanes[i].point());
    }
    r3d_clip(&ico3d, planes, 3);
    r3d_reduce(&ico3d, &vol0, 0);
  }

  // Now clip through our interface and see if we get the same answer.
  const FacetedVolume ico_clip = clipFacetedVolume(ico, clipPlanes);
  if (not fuzzyEqual(ico_clip.volume(), vol0, 1.0e-10)) return "ERROR: clipped volume mismatch: " + to_string(ico_clip.volume()) + " != " + to_string(vol0);
  
  // Must be OK.
  return "OK";
}

}
