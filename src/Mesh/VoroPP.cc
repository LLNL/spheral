//---------------------------------Spheral++----------------------------------//
// VoroPP
//----------------------------------------------------------------------------//
#include <iostream>
#include <sstream>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/tokenizer.hpp"

#include "VoroPP.hh"
#include "MeshConstructionUtilities.hh"
#include "findMatchingVertex.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/boundPointWithinBox.hh"
#include "Utilities/removeElements.hh"

#include "voro++/voro++.cc"

#include "Utilities/timingUtilities.hh"

#define NLAYERS 1

namespace Spheral {
namespace MeshSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Provide a few specialized classes to help reading Voro++ formatted data.
//------------------------------------------------------------------------------
struct VoroVector {
  double x, y, z;
  void fillVector(Dim<2>::Vector& vec) { vec.x(x); vec.y(y); }
  void fillVector(Dim<3>::Vector& vec) { vec.x(x); vec.y(y); vec.z(z); }
};

struct VoroFaceIndices {
  std::vector<unsigned> indices;
};

inline
std::istream&
operator>>(std::istream& is, VoroVector& val) {
  char c;
  is >> c;
  CHECK(c == '(');
  std::string dummy = "";
  double* valptr = &val.x;
  while (c != ')') {
    is >> c;
    if (c == ',' or c == ')') {
      *valptr = boost::lexical_cast<double>(dummy);
      dummy = "";
      ++valptr;
    } else {
      dummy += c;
    }
  }
  return is;
}

inline
std::istream&
operator>>(std::istream& is, VoroFaceIndices& val) {
  val.indices = std::vector<unsigned>();
  char c;
  is >> c;
  CHECK2(c == '(', "c is |" << c << "|");
  std::string dummy = "";
  while (c != ')') {
    is >> c;
    if (c == ',' or c == ')') {
      val.indices.push_back(boost::lexical_cast<unsigned>(dummy));
      dummy = "";
    } else {
      dummy += c;
    }
  }
  return is;
}

//------------------------------------------------------------------------------
// A helper to apply mapping of indices.
//------------------------------------------------------------------------------
inline
void
mapIndex(unsigned& i, const vector<int>& mapIndex, const unsigned maxValue) {
  REQUIRE(i < mapIndex.size());
  REQUIRE(int(mapIndex[i]) < int(maxValue));
  if (mapIndex[i] >= 0) i = mapIndex[i];
}

//------------------------------------------------------------------------------
// VorPP::VoroPP(...)  2-D
//------------------------------------------------------------------------------
template<>
VoroPP<Dim<2> >::
VoroPP(const vector<Dim<2>::Vector>& generators,
       const Dim<2>::Vector& xmin,
       const Dim<2>::Vector& xmax,
       const unsigned nx,
       const unsigned ny,
       const unsigned nz,
       const double edgeTol):
  mNumGenerators(generators.size()),
  mNx(nx),
  mNy(ny),
  mNz(nz),
  mEdgeTol(edgeTol),
  mScale(computeScale(xmin, xmax)),
  mGeneratorsPtr(&generators),
  mXmin(xmin),
  mXmax(xmax),
  mContainerPtr(new container(0.0, (xmax - xmin).x()/mScale,
                              0.0, (xmax - xmin).y()/mScale,
                              0.0, 1.0,
                              mNx, mNy, mNz,       // The number of sub-regions in each dimension.
                              false, false, false, // Periodic?
                              8)) {
  unsigned igen, ilayer, offset;
  const double dz = 1.0/NLAYERS;
  double z;
  Vector gi;

  // Add the generators to the container.
  for (ilayer = 0; ilayer != NLAYERS; ++ilayer) {
    for (igen = 0; igen != mNumGenerators; ++igen) {
      gi = (generators[igen] - xmin) / mScale;
      z = (ilayer + 0.5)*dz;
      offset = ilayer * mNumGenerators;
      mContainerPtr->put(igen + offset, gi.x(), gi.y(), z);
    }
  }
}

//------------------------------------------------------------------------------
// VorPP::VoroPP(...)  3-D
//------------------------------------------------------------------------------
template<>
VoroPP<Dim<3> >::
VoroPP(const vector<Dim<3>::Vector>& generators,
       const Dim<3>::Vector& xmin,
       const Dim<3>::Vector& xmax,
       const unsigned nx,
       const unsigned ny,
       const unsigned nz,
       const double edgeTol):
  mNumGenerators(generators.size()),
  mGeneratorsPtr(&generators),
  mNx(nx),
  mNy(ny),
  mNz(nz),
  mEdgeTol(edgeTol),
  mScale(computeScale(xmin, xmax)),
  mXmin(xmin),
  mXmax(xmax),
  mContainerPtr(new container(0.0, (xmax - xmin).x()/mScale,
                              0.0, (xmax - xmin).y()/mScale,
                              0.0, (xmax - xmin).z()/mScale,
                              mNx, mNy, mNz,       // The number of sub-regions in each dimension.
                              false, false, false, // Periodic?
                              8)) {

  unsigned i, j, k, ijk, igen, ncells = nx*ny*nz;
  Vector gi;

  // Add the generators to the container.
  for (igen = 0; igen != mNumGenerators; ++igen) {
    gi = (generators[igen] - xmin) / mScale;
    subRegion(gi, i, j, k);
    ijk = i + nx*(j + ny*k);
    CHECK(ijk < ncells);
    mContainerPtr->put(igen, gi.x(), gi.y(), gi.z());
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoroPP<Dimension>::
~VoroPP() {
  ENSURE(mContainerPtr.use_count() == 1);
  mContainerPtr->clear();
}

//------------------------------------------------------------------------------
// Return the scaling factor.
//------------------------------------------------------------------------------
template<typename Dimension>
double
VoroPP<Dimension>::
scale() const {
  return mScale;
}

//------------------------------------------------------------------------------
// Add a wall.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoroPP<Dimension>::
addBoundary(MeshWall<Dimension>& meshWall) {
  typedef vector<boost::shared_ptr<wall> > vector_of_walls;
  vector_of_walls walls = meshWall.wallPtrs();
  for (typename vector_of_walls::iterator itr = walls.begin();
       itr != walls.end();
       ++itr) {
    mContainerPtr->add_wall(**itr);
  }
}

//------------------------------------------------------------------------------
// Return the vertex and face description of all cells.
// We use our local Cell class to encapsulate the cell by cell info.
//------------------------------------------------------------------------------
template<>
vector<unsigned>
VoroPP<Dim<2> >::
allCells(vector<Cell<Dim<2> > >& cells) const {
  typedef Dim<2> Dimension;
  const unsigned UNSET = numeric_limits<unsigned>::max();

  // Clear the result.
  cells = vector<Cell<Dimension> >(mNumGenerators);

  // Read out the entire string from Voro++.
  Timing::Time t0 = Timing::currentTime();
  ostringstream oss;
  mContainerPtr->print_all_custom("%i %w %P %s %t %f %C %n", oss);
  const string everything = oss.str();
//   cerr << "Everything:  " << endl
//        << everything << endl;
  if (Process::getRank() == 0) cerr << "VoroPP:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to call print_all_custom with "
                                    << mNumGenerators << " generators." << endl;
  mContainerPtr->clear();

  const double Vscale = mScale*mScale;
  vector<unsigned> flagsMe(mNumGenerators, 0), result;

  // Read out the cell by cell stuff from the string.
  istringstream iss(everything);

  int otherGen;
  unsigned nv, nf, i, j, k, n, igen, numDone = 0U;
  double xc, yc, zc;

  unsigned iigen, ilayer, iBotFace, i0, i1;
  double z0layer, area, volume;
  Vector centroid;
  vector<unsigned> voro2real, realFaceOrder, neighbors;
  vector<vector<unsigned> > faceVertices;
  vector<Vector> vertices;
  vector<VoroVector> voroVerts;
  vector<VoroFaceIndices> voroFaceIndices;
  
  const unsigned maxline = 65536;
  char line[maxline];

  while (iss >> iigen and numDone < mNumGenerators) {
    igen = iigen % mNumGenerators;
    if (flagsMe[igen] == 1) {

      // We've already hit this generator in a different layer, so skip it.
      iss.getline(line, maxline);

    } else {

      CHECK2(igen < mNumGenerators, igen << " " << mNumGenerators);
      CHECK(flagsMe[igen] == 0);
      flagsMe[igen] = 1;
      ++numDone;

      // What layer are we in?
      ilayer = iigen / mNumGenerators;
      CHECK(ilayer < NLAYERS);
      z0layer = ilayer*1.0/NLAYERS;

      // Read the number of vertices.
      iss >> nv;

      // Read out this cell's vertex coordinates.
      CHECK(nv >= 6);
      voroVerts = vector<VoroVector>(nv);
      for (i = 0; i != nv; ++i) {
        iss >> voroVerts[i];
      }
      CHECK(voroVerts.size() == nv);

      // Read out the full (3D) face vertex info.
      iss >> nf;
      CHECK(nf >= 5);
      iBotFace = UNSET;
      voroFaceIndices = vector<VoroFaceIndices>(nf);
      for (i = 0; i != nf; ++i) {
        iss >> voroFaceIndices[i];

        // Identify the face with z = z0layer -- this is our target.
        zc = z0layer;
        j = 0;
        while (j != voroFaceIndices[i].indices.size() and
               abs(voroVerts[voroFaceIndices[i].indices[j]].z - z0layer) < 1.0e-4) ++j;
        if (j == voroFaceIndices[i].indices.size()) {
          CHECK(iBotFace == UNSET);
          iBotFace = i;
        }
      }
      CHECK(iBotFace != UNSET);

      // Pick out the nodes of the bottom face.
      j = 0;
      n = voroFaceIndices[iBotFace].indices.size();
      vertices = vector<Vector>(n);
      faceVertices = vector<vector<unsigned> >(n, vector<unsigned>(2));
      voro2real = vector<unsigned>(nv, UNSET);
      for (i = 0; i != n; ++i) {
        k = voroFaceIndices[iBotFace].indices[i];
        voroVerts[k].fillVector(vertices[i]);
        vertices[i] = vertices[i]*mScale + mXmin;
        faceVertices[i][0] = i;
        faceVertices[i][1] = (i + 1) % n;
        voro2real[k] = i;
      }

      // Extract our real face ordering.
      realFaceOrder = vector<unsigned>(nf, UNSET);
      for (i = 0; i != nf; ++i) {
        vector<unsigned> thpt;
        for (j = 0; j != voroFaceIndices[i].indices.size(); ++j) {
          CHECK(voroFaceIndices[i].indices[j] < nv);
          k = voro2real[voroFaceIndices[i].indices[j]];
          if (k != UNSET) thpt.push_back(k);
        }
        if (thpt.size() == 2) {
          if (thpt[0] < thpt[1]) {
            i0 = thpt[0];
            i1 = thpt[1];
          } else {
            i0 = thpt[1];
            i1 = thpt[0];
          }
          CHECK((i1 == i0 + 1) or (i0 == 0 and i1 == n - 1));
          if (i0 == 0 and i1 == 1) {
            realFaceOrder[i] = 0;
          } else if (i0 == 0 and i1 == n - 1) {
            realFaceOrder[i] = n - 1;
          } else {
            realFaceOrder[i] = i0;
          }
        }
      }

      // Get the face areas (resulting in the zone "volume" in this 2D case).
      for (i = 0; i != nf; ++i) {
        iss >> area;
        if (i == iBotFace) volume = area*Vscale;
      }

      // Get the centroid.
      iss >> xc >> yc >> zc;
      centroid = Vector(xc, yc)*mScale + mXmin;

      // Finally the set of neighbor generators.
      neighbors = vector<unsigned>(n, UNSET);
      for (i = 0; i != nf; ++i) {
        iss >> otherGen;
        if (realFaceOrder[i] != UNSET) {
          CHECK(realFaceOrder[i] < neighbors.size());
          neighbors[realFaceOrder[i]] = (otherGen < 0 ? UNSET : otherGen);
        }
      }

      // Add the Cell to the result.
      cells[igen] = Cell<Dimension>(igen, volume, centroid, vertices, faceVertices, neighbors, mEdgeTol);
    }
  }

  // Check if all generators were created.  If not, return failure.
  if (numDone != mNumGenerators) {
    for (igen = 0; igen != mNumGenerators; ++igen) {
      if (flagsMe[igen] == 0) {
        result.push_back(igen);
        // cerr << "Missed " << igen << endl;
      }
    }
    return result;
  }

  // Otherwise, success!
  return vector<unsigned>();
}

//------------------------------------------------------------------------------
// Return the vertex and face description of all cells.
// We use our local Cell class to encapsulate the cell by cell info.
//------------------------------------------------------------------------------
template<>
vector<unsigned>
VoroPP<Dim<3> >::
allCells(vector<Cell<Dim<3> > >& cells) const {
  typedef Dim<3> Dimension;
  const unsigned UNSET = numeric_limits<unsigned>::max();

  // Clear the result.
  cells = vector<Cell<Dimension> >(mNumGenerators);
  
  // Read out the entire string from Voro++.
  Timing::Time t0 = Timing::currentTime();
  ostringstream oss;
  mContainerPtr->print_all_custom("%i %w %P %s %t %v %C %n", oss);
  const string everything = oss.str();
//   cerr << "Everything:  " << endl
//        << everything << endl;
  if (Process::getRank() == 0) cerr << "VoroPP:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to call print_all_custom with "
                                    << mNumGenerators << " generators." << endl;
  mContainerPtr->clear();

  const double Vscale = mScale*mScale*mScale;
  vector<unsigned> flagsMe(mNumGenerators, 0), result;

  // Read out the cell by cell stuff from the string.
  istringstream iss(everything);

  double volume;
  Vector centroid;
  vector<Vector> vertices;
  vector<vector<unsigned> > faceVertices;
  vector<unsigned> neighbors;

  int otherGen;
  unsigned nv, nf, i, j, k, n, igen, numDone = 0U;
  double xc, yc, zc;
  VoroVector vvec;
  VoroFaceIndices f;

  while (iss >> igen) {
    CHECK2(igen < mNumGenerators, igen << " " << mNumGenerators);
    CHECK(flagsMe[igen] == 0);
    flagsMe[igen] = 1;
    ++numDone;

    // Read the number of vertices.
    iss >> nv;

    // Read out this cell's vertex coordinates.
    CHECK(nv >= 4);
    vertices = vector<Vector>();
    for (i = 0; i != nv; ++i) {
      iss >> vvec;
      vertices.push_back(Vector());
      vvec.fillVector(vertices.back());
      vertices.back() = vertices.back()*mScale + mXmin;
    }
    CHECK(vertices.size() == nv);

    // Read out the face vertex info.
    iss >> nf;
    CHECK(nf >= 4);
    faceVertices = vector<vector<unsigned> >();
    for (i = 0; i != nf; ++i) {
      iss >> f;
      faceVertices.push_back(f.indices);
    }
    CHECK(faceVertices.size() == nf);

    // Get the volume of the cell.
    iss >> volume;
    volume *= Vscale;

    // Get the centroid.
    iss >> xc >> yc >> zc;
    centroid = Vector(xc, yc, zc)*mScale + mXmin;

    // Finally the set of neighbor generators.
    neighbors = vector<unsigned>();
    neighbors.reserve(nf);
    for (i = 0; i != nf; ++i) {
      iss >> otherGen;
      neighbors.push_back(otherGen >= 0 ? otherGen : UNSET);
    }
    CHECK(neighbors.size() == nf);

    // Add the Cell to the result.
    cells[igen] = Cell<Dimension>(igen, volume, centroid, vertices, faceVertices, neighbors, mEdgeTol);
  }

  // Check if all generators were created.  If not, return failure.
  if (numDone != mNumGenerators) {
    for (igen = 0; igen != mNumGenerators; ++igen) {
      if (flagsMe[igen] == 0) {
        result.push_back(igen);
        // cerr << "Missed " << igen << endl;
      }
    }
    return result;
  }

  // Otherwise, success!
  return vector<unsigned>();
}

//------------------------------------------------------------------------------
// Find the indices describing the sub-region containing point p.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoroPP<Dimension>::
subRegion(const typename Dimension::Vector& p, unsigned& i, unsigned& j, unsigned& k) const {
  REQUIRE2(testPointInBox(p, Vector(0,0,0), Vector(1,1,1)),
           "Point not in box:  " << p << " not in [" << mXmin << ", " << mXmax << "]");
  const Vector3d deltaInv(mNx, mNy, mNz);
  i = max(0U, min(mNx - 1, unsigned(p.x()*deltaInv.x())));
  j = max(0U, min(mNy - 1, unsigned(p.y()*deltaInv.y())));
  k = max(0U, min(mNz - 1, unsigned(p.z()*deltaInv.z())));
}

//------------------------------------------------------------------------------
// Compute the scaling factor.
//------------------------------------------------------------------------------
template<typename Dimension>
double
VoroPP<Dimension>::
computeScale(const typename Dimension::Vector& xmin,
             const typename Dimension::Vector& xmax) const {
  const double f = (xmax - xmin).maxElement();
  REQUIRE(f > 1.0e-10);
  return f;
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class VoroPP<Dim<2> >;
template class VoroPP<Dim<3> >;

}
}

// This is *heinous*, but due to the fact that Voro++ will multiply define things
// if you try to include it's headers in different files, we have to build the MeshWall compiled
// bits here in the same .o.
#include "MeshWall.cc"
