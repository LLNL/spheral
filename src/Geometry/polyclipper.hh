//---------------------------------PolyClipper--------------------------------//
// Clip a faceted volume (polygon or polyhedron) by a set of planes in place.
//
// We use the convention that any portion of the faceted volume "below" the 
// plane is clipped, i.e., only the portion of the faceted volume "above" the 
// plane
//    plane.compare(point) >= 0
// is retained.
//
// The algorithms herein are based on R3D as outlined in 
// Powell, D., & Abel, T. (2015). An exact general remeshing scheme applied to 
// physically conservative voxelization. Journal of Computational Physics, 297, 340–356.
//
// Created by J. Michael Owen, Tue Nov 28 10:00:51 PST 2017
//----------------------------------------------------------------------------//
#ifndef __PolyClipper_hh__
#define __PolyClipper_hh__

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

#include <string>
#include <vector>
#include <limits>

namespace PolyClipper {

//------------------------------------------------------------------------------
// 2D plane.
//------------------------------------------------------------------------------
struct Plane2d {
  typedef Spheral::Dim<2>::Vector Vector;
  double dist;                       // Signed distance to the origin
  Vector normal;                     // Unit normal
  int ID;                            // ID for the plane, used to label vertices
  Plane2d()                                                  : dist(0.0), normal(1,0), ID(std::numeric_limits<int>::min()) {}
  Plane2d(const double d, const Vector& nhat)                : dist(d), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane2d(const Vector& p, const Vector& nhat)               : dist(-p.dot(nhat)), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane2d(const Vector& p, const Vector& nhat, const int id) : dist(-p.dot(nhat)), normal(nhat), ID(id) {}
  Plane2d(const Plane2d& rhs)                                : dist(rhs.dist), normal(rhs.normal), ID(rhs.ID) {}
  Plane2d& operator=(const Plane2d& rhs)                     { dist = rhs.dist; normal = rhs.normal; ID = rhs.ID; return *this; }
  bool operator==(const Plane2d& rhs) const                  { return (dist == rhs.dist and normal == rhs.normal); }
  bool operator!=(const Plane2d& rhs) const                  { return not (*this == rhs); }
  bool operator< (const Plane2d& rhs) const                  { return (dist < rhs.dist); }
  bool operator> (const Plane2d& rhs) const                  { return (dist > rhs.dist); }
};

//------------------------------------------------------------------------------
// 3D plane.
//------------------------------------------------------------------------------
struct Plane3d {
  typedef Spheral::Dim<3>::Vector Vector;
  Vector normal;                     // Unit normal
  double dist;                       // Signed distance to the origin
  int ID;                            // ID for the plane, used to label vertices
  Plane3d()                                                  : dist(0.0), normal(1,0,0), ID(std::numeric_limits<int>::min()) {}
  Plane3d(const double d, const Vector& nhat)                : dist(d), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane3d(const Vector& p, const Vector& nhat)               : dist(-p.dot(nhat)), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane3d(const Vector& p, const Vector& nhat, const int id) : dist(-p.dot(nhat)), normal(nhat), ID(id) {}
  Plane3d(const Plane3d& rhs)                                : dist(rhs.dist), normal(rhs.normal), ID(rhs.ID) {}
  Plane3d& operator=(const Plane3d& rhs)                     { dist = rhs.dist; normal = rhs.normal; ID = rhs.ID; return *this; }
  bool operator==(const Plane3d& rhs) const                  { return (dist == rhs.dist and normal == rhs.normal); }
  bool operator!=(const Plane3d& rhs) const                  { return not (*this == rhs); }
  bool operator< (const Plane3d& rhs) const                  { return (dist < rhs.dist); }
  bool operator> (const Plane3d& rhs) const                  { return (dist > rhs.dist); }
};

//------------------------------------------------------------------------------
// The 2D vertex struct, which we use to encode polygons.
//------------------------------------------------------------------------------
struct Vertex2d {
  typedef Spheral::Dim<2>::Vector Vector;
  Vector position;
  std::pair<int, int> neighbors;
  int comp;
  mutable int ID;                            // convenient, but sneaky
  mutable std::set<int> clips;               // the planes (if any) that created this point
  Vertex2d():                                  position(),    neighbors(), comp(1), ID(-1), clips() {}
  Vertex2d(const Vector& pos):                 position(pos), neighbors(), comp(1), ID(-1), clips() {}
  Vertex2d(const Vector& pos, const int c)   : position(pos), neighbors(), comp(c), ID(-1), clips() {}
  Vertex2d(const Vertex2d& rhs)              : position(rhs.position), neighbors(rhs.neighbors), comp(rhs.comp), ID(rhs.ID), clips(rhs.clips) {}
  Vertex2d& operator=(const Vertex2d& rhs)   { position = rhs.position; neighbors = rhs.neighbors; comp = rhs.comp; ID = rhs.ID; clips = rhs.clips; return *this; }
  bool operator==(const Vertex2d& rhs) const {
    return (position  == rhs.position and
            neighbors == rhs.neighbors and
            comp      == rhs.comp and
            ID        == rhs.ID);
  }
};

//------------------------------------------------------------------------------
// The 3D vertex struct, which we use to encode polyhedra.
//------------------------------------------------------------------------------
struct Vertex3d {
  typedef Spheral::Dim<3>::Vector Vector;
  Vector position;
  std::vector<int> neighbors;
  int comp;
  mutable int ID;                            // convenient, but sneaky
  mutable std::set<int> clips;               // the planes (if any) that created this point
  Vertex3d():                                  position(),    neighbors(), comp(1), ID(-1), clips() {}
  Vertex3d(const Vector& pos):                 position(pos), neighbors(), comp(1), ID(-1), clips() {}
  Vertex3d(const Vector& pos, const int c):    position(pos), neighbors(), comp(c), ID(-1), clips() {}
  Vertex3d(const Vertex3d& rhs)              : position(rhs.position), neighbors(rhs.neighbors), comp(rhs.comp), ID(rhs.ID), clips(rhs.clips) {}
  Vertex3d& operator=(const Vertex3d& rhs)   { position = rhs.position; neighbors = rhs.neighbors; comp = rhs.comp; ID = rhs.ID; clips = rhs.clips; return *this; }
  bool operator==(const Vertex3d& rhs) const {
    return (position  == rhs.position and
            neighbors == rhs.neighbors and
            comp      == rhs.comp and
            ID        == rhs.ID);
  }
};

//------------------------------------------------------------------------------
// 2D (polygon) methods.
//------------------------------------------------------------------------------
typedef std::vector<Vertex2d> Polygon;

void initializePolygon(Polygon& poly,
                       const std::vector<Spheral::Dim<2>::Vector>& positions,
                       const std::vector<std::vector<int>>& neighbors);

std::string polygon2string(const Polygon& poly);

void convertToPolygon(Polygon& polygon,
                      const Spheral::Dim<2>::FacetedVolume& Spheral_polygon);

std::vector<std::set<int>> convertFromPolygon(Spheral::Dim<2>::FacetedVolume& Spheral_polygon,
                                              const Polygon& polygon);

void moments(double& zerothMoment, Spheral::Dim<2>::Vector& firstMoment,
             const Polygon& polygon);

void clipPolygon(Polygon& poly,
                 const std::vector<Plane2d>& planes);

void collapseDegenerates(Polygon& poly,
                         const double tol);

  std::vector<std::vector<int>> splitIntoTriangles(const Polygon& poly,
                                                   const double tol = 0.0);

//------------------------------------------------------------------------------
// 3D (polyhedron) methods.
//------------------------------------------------------------------------------
typedef std::vector<Vertex3d> Polyhedron;

std::vector<std::vector<int>> extractFaces(const Polyhedron& poly);

void initializePolyhedron(Polyhedron& poly,
                          const std::vector<Spheral::Dim<3>::Vector>& positions,
                          const std::vector<std::vector<int>>& neighbors);

std::string polyhedron2string(const Polyhedron& poly);

void convertToPolyhedron(Polyhedron& polyhedron,
                         const Spheral::Dim<3>::FacetedVolume& Spheral_polyhedron);

std::vector<std::set<int>>  convertFromPolyhedron(Spheral::Dim<3>::FacetedVolume& Spheral_polyhedron,
                                                  const Polyhedron& polyhedron);

void moments(double& zerothMoment, Spheral::Dim<3>::Vector& firstMoment,
             const Polyhedron& polyhedron);

void clipPolyhedron(Polyhedron& poly,
                    const std::vector<Plane3d>& planes);

void collapseDegenerates(Polyhedron& poly,
                         const double tol);

std::vector<std::vector<int>> splitIntoTetrahedra(const Polyhedron& poly, 
                                                  const double tol = 0.0);
}

#endif

