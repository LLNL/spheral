//---------------------------------Spheral++----------------------------------//
// VoroPP2d
// 
// An interface/container class to encapsulate our use of the Voro++ library.
// This one is specialized for the new 2D library.
// This is necessary 'cause Voro++ doesn't build itself, you're supposed to 
// include it's .cc file into your own file.  I can't do this multiple times
// everywhere I use Voro++ because that will cause mutiple instantiations, so 
// hide everything behind a single instantiation in this object.
//
// Created by JMO, Thu Jul 21 11:14:11 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_VoroPP2d__
#define __Spheral_VoroPP2d__

#include <vector>
#include "boost/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

// Forward declare voro++ stuff.
class container_2d;
class voronoicell_2d;

namespace Spheral {
namespace MeshSpace {

//------------------------------------------------------------------------------
// VoroPP2d -- encapsulate the Voro++ container interface.
//------------------------------------------------------------------------------
class VoroPP2d {
  //--------------------------- Public Interface ---------------------------//
public:
  typedef Dim<2> Dimension;
  typedef Dimension::Vector Vector;
  typedef Dimension::FacetedVolume FacetedVolume;
  typedef GeomPlane<Dimension> Plane;

  // Constructors, destructor.
  VoroPP2d(const std::vector<Vector>& generators, const Vector& xmin, const Vector& xmax,
           const unsigned nx = 20,
           const unsigned ny = 20);
  ~VoroPP2d();

  // Return the vertices for all cells.
  // Vertex positions are listed in counter-clockwise order for each cell.
  // Note we assign a unique ID to each vertex, but the vertices themselves are
  // degenerate!
  void allCells(std::vector<Vector>& vertices,
                std::vector<std::vector<unsigned> >& cellVertexIndices) const;

  // Compute which sub-region the given position is in.
  void subRegion(const Vector& p, unsigned& i, unsigned& j) const;

private:
  unsigned mNumGenerators, mNx, mNy;
  const std::vector<Vector>* mGeneratorsPtr;
  Vector mXmin, mXmax;
  boost::shared_ptr<container_2d> mContainerPtr;
  std::vector<unsigned> mGen2CellGen, mVoro2GenOrder;

  // Internal method for common constructor operations.
  void construct(const std::vector<Vector>& generators,
                 const unsigned nx,
                 const unsigned ny);

  // Forbidden methods.
  VoroPP2d();
  VoroPP2d(const VoroPP2d&);
  VoroPP2d& operator=(const VoroPP2d&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace MeshSpace {
    class VoroPP2d;
  }
}

#endif
