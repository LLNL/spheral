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
using namespace boost;

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
  mGeneratorsPtr(&generators),
  mXmin(xmin),
  mXmax(xmax),
  mContainerPtr(new container(xmin.x(), xmax.x(),
                              xmin.y(), xmax.y(),
                              0.0, 1.0, // (xmax.x() - xmin.x())/mNx,
                              mNx, mNy, mNz,       // The number of sub-regions in each dimension.
                              false, false, false, // Periodic?
                              8)) {
  unsigned igen, ilayer, offset;
  const double dz = 1.0/NLAYERS;
  double z;

//   const Vector box = xmax - xmin;
//   const Vector boxInv = Vector(1.0/box.x(),
//                                1.0/box.y());

  // Add the generators to the container.
  for (ilayer = 0; ilayer != NLAYERS; ++ilayer) {
    z = (ilayer + 0.5)*dz;
    offset = ilayer * mNumGenerators;
    for (igen = 0; igen != mNumGenerators; ++igen) {
      mContainerPtr->put(igen + offset, generators[igen].x(), generators[igen].y(), z);
//                          (generators[igen].x() - mXmin.x())*boxInv.x(),
//                          (generators[igen].y() - mXmin.y())*boxInv.y(),
//                          z);
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
  mXmin(xmin),
  mXmax(xmax),
  mContainerPtr(new container(xmin.x(), xmax.x(),
                              xmin.y(), xmax.y(),
                              xmin.z(), xmax.z(),
                              mNx, mNy, mNz,       // The number of sub-regions in each dimension.
                              false, false, false, // Periodic?
                              8)) {

  unsigned i, j, k, ijk, igen, ncells = nx*ny*nz;

  // Add the generators to the container.
  for (igen = 0; igen != mNumGenerators; ++igen) {
    subRegion(generators[igen], i, j, k);
    ijk = i + nx*(j + ny*k);
    CHECK(ijk < ncells);
    mContainerPtr->put(igen, generators[igen].x(), generators[igen].y(), generators[igen].z());
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
// Note we assign a unique ID to each vertex, but the vertices themselves are
// degenerate!
//------------------------------------------------------------------------------
template<>
vector<unsigned>
VoroPP<Dim<2> >::
allCells(vector<vector<Dim<2>::Vector> >& cellVertices,
         vector<vector<vector<unsigned> > >& cellFaceVertices,
         vector<double>& cellVolumes,
         vector<Dim<2>::Vector>& cellCentroids,
         vector<vector<unsigned> >& neighborCells) const {

  // Size the results appropriately.
  cellVertices = vector<vector<Vector> >(mNumGenerators);
  cellFaceVertices = vector<vector<vector<unsigned> > >(mNumGenerators);
  cellVolumes = vector<double>(mNumGenerators);
  cellCentroids = vector<Vector>(mNumGenerators);
  neighborCells = vector<vector<unsigned> >(mNumGenerators);
  
  const unsigned UNSET = numeric_limits<unsigned>::max();
//   const double zmin = -(mXmax.x() - mXmin.x())/mNx;

  // Read out the entire string from Voro++.
//   Timing::Time t0 = Timing::currentTime();
  ostringstream oss;
  mContainerPtr->print_all_custom("%i %w %P %s %t %f %C %n", oss);
  const string everything = oss.str();
//   if (Process::getRank() == 0) cerr << "VoroPP:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to call print_all_custom with "
//                                     << mNumGenerators << " generators." << endl;
  mContainerPtr->clear();

  vector<unsigned> flagsMe(mNumGenerators, 0), result;

  // Read out the cell by cell stuff from the string.
  istringstream iss(everything);
  const unsigned maxline = 65536;
  char line[maxline];

  bool useGen, killEdge, done;
  int otherGen;
  unsigned nv, nf, i, j, k, n, igen, jgen, iigen, ilayer, iBotFace, numDone = 0U;
  double xc, yc, zc, area, maxEdge, z0layer;
  vector<unsigned> voro2real, edges2kill;
  vector<int> realFaceOrder, newVertexOrder;
  vector<VoroVector> voroVerts;
  vector<VoroFaceIndices> voroFaceIndices;
  VoroVector vvec;
  VoroFaceIndices f;
  while (iss >> iigen and numDone < mNumGenerators) {
    igen = iigen % mNumGenerators;
    if (flagsMe[igen] == 1) {
      // We've already hit this generator in a different layer, so skip it.
      iss.getline(line, maxline);

    } else {
      CHECK2(igen < mNumGenerators, igen << " " << mNumGenerators);
      flagsMe[igen] = 1;
      ++numDone;
      ilayer = iigen / mNumGenerators;
      CHECK(ilayer < NLAYERS);
      z0layer = ilayer*1.0/NLAYERS;

      // Read the number of vertices.
      iss >> nv;

      // Read out this cell's Voro++ vertex coordinates.
      CHECK(nv >= 6);
      voroVerts = vector<VoroVector>();
      for (i = 0; i != nv; ++i) {
        voroVerts.push_back(VoroVector());
        iss >> voroVerts.back();
      }
      CHECK(voroVerts.size() == nv);

      // Read out the full (3D) face vertex info.
      iss >> nf;
      CHECK(nf >= 5);
      iBotFace = UNSET;
      voroFaceIndices = vector<VoroFaceIndices>();
      for (i = 0; i != nf; ++i) {
        voroFaceIndices.push_back(VoroFaceIndices());
        iss >> voroFaceIndices.back();

        // Identify the face with all zero z coordinates.
        // This represents our 2D cell.
        zc = z0layer;
        for (j = 0; j != voroFaceIndices.back().indices.size(); ++j) {
          zc = max(zc, voroVerts[voroFaceIndices.back().indices[j]].z);
        }
        if (zc < z0layer + 0.01) {
          CHECK(iBotFace == UNSET);
          iBotFace = i;
        }
      }
      CHECK(iBotFace != UNSET);

      // Pick out the nodes of the bottom face.
      j = 0;
      n = voroFaceIndices[iBotFace].indices.size();
      cellVertices[igen].resize(n);
      voro2real = vector<unsigned>(nv, UNSET);
      for (i = 0; i != n; ++i) {
        k = voroFaceIndices[iBotFace].indices[i];
        voroVerts[k].fillVector(cellVertices[igen][i]);
        boundPointWithinBox(cellVertices[igen][i], Vector(0.0,0.0), Vector(1.0, 1.0));
        voro2real[k] = i;
      }

      // Extract our real face ordering.
      realFaceOrder = vector<int>(nf, -1);
      for (i = 0; i != nf; ++i) {
        vector<unsigned> thpt;
        for (j = 0; j != voroFaceIndices[i].indices.size(); ++j) {
          CHECK(voroFaceIndices[i].indices[j] < nv);
          k = voro2real[voroFaceIndices[i].indices[j]];
          if (k != UNSET) thpt.push_back(k);
        }
        if (thpt.size() == 2) {
          if (min(thpt[0], thpt[1]) == 0) {
            if (max(thpt[0], thpt[1]) == 1) {
              realFaceOrder[i] = 0;
            } else {
              realFaceOrder[i] = max(thpt[0], thpt[1]);
            }
          } else {
            realFaceOrder[i] = min(thpt[0], thpt[1]);
          }
        }
      }
      n = cellVertices[igen].size();
      for (i = 0; i != n; ++i) {
        vector<unsigned> thpt;
        thpt.push_back(i);
        thpt.push_back((i + 1) % n);
        cellFaceVertices[igen].push_back(thpt);
      }
      CHECK(cellFaceVertices[igen].size() >= 3);
      // if (not (cellFaceVertices[igen].size() == cellVertices[igen].size())) {
      //   cerr << "Vertices: ";
      //   for (i = 0; i != cellVertices[igen].size(); ++i) cerr << cellVertices[igen][i] << " ";
      //   cerr << endl
      //        << "          ";
      //   for (i = 0; i != voroFaceIndices[iBotFace].indices.size(); ++i) {
      //     k = voroFaceIndices[iBotFace].indices[i];
      //     cerr << "(" << voroVerts[k].x << " " << voroVerts[k].y << " " << voroVerts[k].z << ") ";
      //   }
      //   cerr << endl
      //        << "          ";
      //   for (i = 0; i != nv; ++i) {
      //     cerr << "(" << voroVerts[k].x << " " << voroVerts[k].y << " " << voroVerts[k].z << ") ";
      //   }
      //   cerr << endl
      //        << "          ";
      //   for (i = 0; i != nv; ++i) cerr << voro2real[i] << " ";
      //   cerr << endl
      //        << "Face vertices: ";
      //   for (i = 0; i != cellFaceVertices[igen].size(); ++i) cerr << "(" << cellFaceVertices[igen][i][0] << " " << cellFaceVertices[igen][i][1] << ") ";
      //   cerr << endl;
      // }
      CHECK(cellFaceVertices[igen].size() == cellVertices[igen].size());

      // Get the face areas (resulting in the zone "volume" in this 2D case).
      for (i = 0; i != nf; ++i) {
        iss >> area;
        if (i == iBotFace) cellVolumes[igen] = area;
      }

      // Get the centroid.
      iss >> xc >> yc >> zc;
      cellCentroids[igen] = Vector(xc, yc);

      // Finally the set of neighbor generators.
      neighborCells[igen] = vector<unsigned>(n);
      for (i = 0; i != nf; ++i) {
        iss >> otherGen;
        if (realFaceOrder[i] >= 0) {
          CHECK(realFaceOrder[i] < n);
          neighborCells[igen][realFaceOrder[i]] = (otherGen >= 0 ? (otherGen % mNumGenerators) : UNSET);
        }
      }
    }
  }

  // Check if all generators were created.  If not, return failure.
  if (numDone != mNumGenerators) {
    for (igen = 0; igen != mNumGenerators; ++igen) {
      if (flagsMe[igen] == 0) result.push_back(igen);
    }
    return result;
  }

  // Now do some pruning.
  unsigned ji, jj, nj;
  vector<unsigned>::const_iterator killItr;
  vector<vector<Vector> > newCellVertices;
  done = false;
  while (not done) {
    done = true;
    vector<vector<unsigned> > edges2kill(mNumGenerators);
    for (igen = 0; igen != mNumGenerators; ++igen) {
      // Scan the edges looking for any that fall below the threshold for removal.
      n = cellVertices[igen].size();
      maxEdge = 0.0;
      for (i = 0; i != n; ++i) {
        j = (i + 1) % n;
        maxEdge = max(maxEdge, (cellVertices[igen][i] - cellVertices[igen][j]).magnitude());
      }
      CHECK2(maxEdge > 0.0, igen << " " << maxEdge << " " << (*mGeneratorsPtr)[igen]);
      for (i = 0; i != n; ++i) {
        jgen = neighborCells[igen][i];
        j = (i + 1) % n;
        killEdge = ((cellVertices[igen][i] - cellVertices[igen][j]).magnitude()/maxEdge < mEdgeTol);
        if ((not killEdge) and (jgen != UNSET)) {
          // Check if the neighbor zone across this face knows about us.  If not, kill the edge.
          killEdge = (count(neighborCells[jgen].begin(), neighborCells[jgen].end(), igen) == 0);
        }
        if (killEdge) {
          edges2kill[igen].push_back(i);
          if (jgen != UNSET) {
            k = distance(neighborCells[jgen].begin(), find(neighborCells[jgen].begin(), neighborCells[jgen].end(), igen));
            if (k < neighborCells[jgen].size()) {
              edges2kill[jgen].push_back(k);
            }
          }
          done = false;
        }
      }
    }

    if (not done) {
      newCellVertices = cellVertices;
      // Make a first pass on any edges we're removing and make sure everyone agrees about
      // new node positions.
      for (igen = 0; igen != mNumGenerators; ++igen) {
        if (edges2kill[igen].size() > 0) {
          newCellVertices[igen] = cellVertices[igen];
          sort(edges2kill[igen].begin(), edges2kill[igen].end());
          edges2kill[igen].erase(unique(edges2kill[igen].begin(), edges2kill[igen].end()), edges2kill[igen].end());
          n = cellVertices[igen].size();
          for (killItr = edges2kill[igen].begin(); killItr != edges2kill[igen].end(); ++killItr) {
            i = *killItr;
            CHECK(i < n);
            j = (i + 1) % n;
            jgen = neighborCells[igen][i];
            
            // Does the neighbor really line up with us?
            if (jgen != UNSET) {
              nj = cellVertices[jgen].size();
              ji = findMatchingVertex(i, cellVertices[igen], cellVertices[jgen]);
              CHECK(ji < nj);
              jj = (ji == 0 ? nj - 1 : ji - 1);
              if ((cellVertices[igen][i] - cellVertices[jgen][ji]).magnitude2() +
                  (cellVertices[igen][j] - cellVertices[jgen][jj]).magnitude2() > 1.0e-8) jgen = UNSET;
            }
            if (jgen == UNSET) {
              newCellVertices[igen][j] = 0.5*(cellVertices[igen][i] + cellVertices[igen][j]);
              newCellVertices[igen][i] = newCellVertices[igen][j];
            } else {
              newCellVertices[igen][j] = 0.25*(cellVertices[igen][i] + cellVertices[igen][j] +
                                               cellVertices[jgen][ji] + cellVertices[jgen][jj]);
              newCellVertices[igen][i] = newCellVertices[igen][j];
              newCellVertices[jgen][ji] == newCellVertices[igen][j];
              newCellVertices[jgen][jj] == newCellVertices[igen][j];
            }
          }
        }
      }

      // Now remove edges & nodes.
      for (igen = 0; igen != mNumGenerators; ++igen) {
        if (edges2kill[igen].size() > 0) {
          n = cellVertices[igen].size();
          vector<int> newVertexOrder(n, -1);
          killItr = edges2kill[igen].begin();
          k = 0;
          for (i = 0; i != n; ++i) {
            if (killItr < edges2kill[igen].end() and i == *killItr) {
              ++killItr;
            } else {
              newVertexOrder[i] = k++;
            }
          }
          CHECK2(killItr == edges2kill[igen].end(), killItr - edges2kill[igen].begin() << " " << edges2kill[igen].size());
          CHECK2(k == (n - edges2kill[igen].size()), k << " " << n << " " << edges2kill[igen].size());
          CHECK2(k >= 3, n << " " << edges2kill[igen].size());
          for (i = 0; i != n; ++i) {
            CHECK(cellFaceVertices[igen][i].size() == 2);
            mapIndex(cellFaceVertices[igen][i][0], newVertexOrder, k);
            mapIndex(cellFaceVertices[igen][i][1], newVertexOrder, k);
          }
          removeElements(newCellVertices[igen], edges2kill[igen]);
          removeElements(cellFaceVertices[igen], edges2kill[igen]);
          removeElements(neighborCells[igen], edges2kill[igen]);
          CHECK(newCellVertices[igen].size() == k);
          CHECK(cellFaceVertices[igen].size() == k);
          CHECK(neighborCells[igen].size() == k);
        }
      }
      cellVertices = newCellVertices;
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
//     // Check that all generators were computed.
//     const unsigned numDone = accumulate(flagsMe.begin(), flagsMe.end(), 0U);
//     if (numDone != mNumGenerators) {
//       cerr << "Box range:  " << mXmin << " " << mXmax << endl;
//       for (igen = 0; igen != NLAYERS * mNumGenerators; ++igen) {
//         if (flagsMe[igen] == 0U) cerr << "  --> " << igen << " " << (*mGeneratorsPtr)[igen % mNumGenerators] << endl;
//       }
//     }
//     ENSURE2(numDone == mNumGenerators, "VoroPP: computed " << numDone << " of " << mNumGenerators << " generators.");
    
    // Look for symmetry in neighbors.
    unsigned jgen;
    for (igen = 0; igen != mNumGenerators; ++igen) {
      nf = cellVertices[igen].size();
//       if (not (cellFaceVertices[igen].size() == nf)) {
//         cerr << "Vertices: ";
//         for (i = 0; i != nf; ++i) cerr << cellVertices[igen][i] << " ";
//         cerr << endl
//              << "          ";
//         for (i = 0; i != voroFaceIndices[iBotFace].indices.size(); ++i) {
//           k = voroFaceIndices[iBotFace].indices[i];
//           cerr << "(" << voroVerts[k].x << " " << voroVerts[k].y << " " << voroVerts[k].z << ") ";
//         }
//         cerr << endl
//              << "Face vertices: ";
//         for (i = 0; i != cellFaceVertices[igen].size(); ++i) cerr << "(" << cellFaceVertices[igen][i][0] << " " << cellFaceVertices[igen][i][1] << ") ";
//         cerr << "Cell centroid: " << "(" << xc << " " << yc << " " << zc << ")" << endl;
//         cerr << endl;
//       }
      ENSURE2(cellFaceVertices[igen].size() == nf, cellFaceVertices[igen].size() << " " << nf);
      ENSURE(neighborCells[igen].size() == nf);
      for (i = 0; i != nf; ++i) {
        jgen = neighborCells[igen][i];
        ENSURE2(jgen == UNSET or count(neighborCells[igen].begin(), neighborCells[igen].end(), jgen) == 1,
                igen << " " << jgen << " : " << count(neighborCells[igen].begin(), neighborCells[igen].end(), jgen));
        if (jgen != UNSET) {
//           if (!(count(neighborCells[jgen].begin(), neighborCells[jgen].end(), igen) == 1)) {
//             cerr << "Blago!  " << endl
//                  << "Generator " << igen << endl
//                  << "    ";
//             for (j = 0; j != cellVertices[igen].size(); ++j) cerr << cellVertices[igen][j] << " ";
//             cerr << endl
//                  << "    ";
//             for (j = 0; j != neighborCells[igen].size(); ++j) cerr << neighborCells[igen][j] << " ";
//             cerr << endl
//                  << "Generator " << jgen << endl
//                  << "    ";
//             for (j = 0; j != cellVertices[jgen].size(); ++j) cerr << cellVertices[jgen][j] << " ";
//             cerr << endl
//                  << "    ";
//             for (j = 0; j != neighborCells[jgen].size(); ++j) cerr << neighborCells[jgen][j] << " ";
//             cerr << endl;
//           }
          ENSURE2(count(neighborCells[jgen].begin(), neighborCells[jgen].end(), igen) == 1,
                  jgen << " " << igen << " : " << count(neighborCells[jgen].begin(), neighborCells[jgen].end(), igen));
        }
        j = (i + 1) % nf;
        ENSURE(cellFaceVertices[igen][i] != cellFaceVertices[igen][j]);
      }
    }
  }
  END_CONTRACT_SCOPE;

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Return the vertex and face description of all cells.
// Note we assign a unique ID to each vertex, but the vertices themselves are
// degenerate!
//------------------------------------------------------------------------------
template<>
vector<unsigned>
VoroPP<Dim<3> >::
allCells(vector<Cell<Dim<3> > >& cells) const {
  typedef Dim<3> Dimension;

  // Clear the result.
  cells = vector<Cell<Dimension> >();
  
  const unsigned UNSET = numeric_limits<unsigned>::max();

  // Read out the entire string from Voro++.
  Timing::Time t0 = Timing::currentTime();
  ostringstream oss;
  mContainerPtr->print_all_custom("%i %w %P %s %t %v %C %n", oss);
  const string everything = oss.str();
  // cerr << "Everything:  " << endl
  //      << everything << endl;
  if (Process::getRank() == 0) cerr << "VoroPP:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to call print_all_custom with "
                                    << mNumGenerators << " generators." << endl;
  mContainerPtr->clear();

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

    // Get the centroid.
    iss >> xc >> yc >> zc;
    centroid = Vector(xc, yc, zc);

    // Finally the set of neighbor generators.
    neighbors = vector<unsigned>();
    neighbors.reserve(nf);
    for (i = 0; i != nf; ++i) {
      iss >> otherGen;
      neighbors.push_back(otherGen >= 0 ? otherGen : UNSET);
    }
    CHECK(neighbors.size() == nf);

    // Add the Cell to the result.
    cells.push_back(Cell<Dimension>(igen, volume, centroid, vertices, faceVertices, neighbors, mEdgeTol));
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
  REQUIRE2(testPointInBox(p, mXmin, mXmax),
           "Point not in box:  " << p << " not in [" << mXmin << ", " << mXmax << "]");
  const Vector3d deltaInv(mNx*safeInv(mXmax.x() - mXmin.x()),
                          mNy*safeInv(mXmax.y() - mXmin.y()),
                          mNz*safeInv(mXmax.z() - mXmin.z()));
  i = max(0U, min(mNx - 1, unsigned((p.x() - mXmin.x())*deltaInv.x())));
  j = max(0U, min(mNy - 1, unsigned((p.y() - mXmin.y())*deltaInv.y())));
  k = max(0U, min(mNz - 1, unsigned((p.z() - mXmin.z())*deltaInv.z())));
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
