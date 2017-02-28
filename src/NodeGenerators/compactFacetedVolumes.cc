#include <algorithm>

#include "compactFacetedVolumes.hh"

#include "NodeList/FluidNodeList.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/GammaLawGas.hh"
#include "DataBase/DataBase.hh"
#include "Utilities/packElement.hh"
#include "Utilities/allReduce.hh"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Push FacetedVolume shapes together inside a surface, but excluding mutual
// overlap.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned compactFacetedVolumes(std::vector<typename Dimension::FacetedVolume>& shapes,
                               std::vector<typename Dimension::Vector>& centers,
                               std::vector<int>& flags,
                               const typename Dimension::FacetedVolume& surface,
                               const double depthmax,
                               const unsigned surfaceIterations,
                               const unsigned maxIterations,
                               const double dispfrac,
                               const double maxoverlapfrac) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  using FieldSpace::Field;

  // Preconditions.
  const unsigned nshapes = shapes.size();
  VERIFY(centers.size() == nshapes);
  VERIFY(flags.size() == nshapes);

  // Only proceed if there's work to do!
  int flagmax = *max_element(flags.begin(), flags.end());
  if (allReduce(flagmax, MPI_MAX, Communicator::communicator()) != 2) return 0;

  // Carve up the shapes range in parallel.
  const size_t rank = Process::getRank();
  const size_t nprocs = Process::getTotalNumberOfProcesses();
  const size_t ndomain0 = nshapes/nprocs;
  const size_t remainder = nshapes % nprocs;
  CHECK(remainder < nprocs);
  size_t ndomain = ndomain0;
  if (rank < remainder) ndomain += 1;
  const size_t imin = rank*ndomain0 + min(rank, remainder);
  const size_t imax = imin + ndomain;

  // Build a fake NodeList with Neighbor to help selecting shapes that might be neighbors.
  const double length = (surface.xmax() - surface.xmin()).maxElement();
  size_t numActiveShapes = 0;
  for (size_t i = 0; i != nshapes; ++i) {
    if (flags[i] != 0) ++numActiveShapes;
  }

  // Make a temporary NodeList so we can use it's neighbor logic.
  Material::PhysicalConstants constants(1.0, 1.0, 1.0);
  Material::GammaLawGas<Dimension> eos(2.0, 2.0, constants, 0.0, 1e10, Material::MaterialPressureMinType::PressureFloor);
  NodeSpace::FluidNodeList<Dimension> nodes("shapes", eos, numActiveShapes, 0,
                                            1e-10, 1e10, 1.0, 1.01, 500, 0.0, 1e10);
  NeighborSpace::NestedGridNeighbor<Dimension> neighbor(nodes, NeighborSpace::NeighborSearchType::GatherScatter, 31, 2.8*length, Vector::zero, 1.0, 1);
  nodes.registerNeighbor(neighbor);
  DataBaseSpace::DataBase<Dimension> db;
  db.appendNodeList(nodes);
  Field<Dimension, Vector>& pos = nodes.positions();
  Field<Dimension, SymTensor>& H = nodes.Hfield();
  Field<Dimension, double> radius("radius", nodes);

  // Figure out the effective size per shape.
  size_t j = 0;
  vector<size_t> shape2node;
  for (size_t i = 0; i < nshapes; ++i) {
    shape2node.push_back(j);
    if (flags[i] != 0) {
      pos[j] = centers[i];
      double Ri = 0.0;
      if (flags[i] == 1) {
        for (const auto& v: shapes[j].vertices()) {
          Ri = std::max(Ri, (v - centers[i]).magnitude());
        }
      } else {
        for (const auto& v: shapes[j].vertices()) {
          Ri = std::max(Ri, v.magnitude());
        }
      }
      radius[j] = Ri;
      H[j] = SymTensor::one / (2.0*Ri);
      ++j;
    }
  }
  CHECK(shape2node.size() == nshapes);

  // First iterate the requested number of iterations with surface repulsion.
  unsigned iter = 0;
  double maxoverlap = 10.0*maxoverlapfrac;
  while (iter < surfaceIterations or (iter < maxIterations and maxoverlap > maxoverlapfrac)) {
    iter += 1;
    vector<Vector> displacements(nshapes);

    // Update the neighbor info.
    neighbor.updateNodes();
    db.updateConnectivityMap(false);
    const NeighborSpace::ConnectivityMap<Dimension>& cm = db.connectivityMap();

    // Add the repulsive component from the asteroid surface.
    for (auto i = imin; i < imax; ++i) {
      if (flags[i] == 2) {
        const Vector p = surface.closestPoint(centers[i]);
        const double d = surface.distance(centers[i]);
        if (surface.contains(centers[i])) {
          displacements[i] += dispfrac*(1.0 - double(iter)/surfaceIterations)*max(0.01, 1.0 - d/depthmax)*Dimension::rootnu(shapes[i].volume()/M_PI)*(centers[i] - p).unitVector();
        } else {
          displacements[i] -= dispfrac*(1.0 - double(iter)/surfaceIterations)*Dimension::rootnu(shapes[i].volume()/M_PI)*(centers[i] - p).unitVector();
        }
      }
    }

    // Look for any overlapping shapes and drive them apart.
    for (auto i = imin; i < imax; ++i) {
      if (flags[i] == 3) {
        const FacetedVolume bi = shapes[i] + centers[i];
        const vector<vector<int>>& neighbors = cm.connectivityForNode(0, shape2node[i]);
        CHECK(neighbors.size() == 1);
        for (auto j: neighbors[0]) {
          CHECK(j != shape2node[i]);
          if ((flags[j] == 1 and bi.intersect(shapes[j])) or
              (flags[j] >= 2 and bi.intersect(shapes[j] + centers[j]))) {
            const double Vi = shapes[i].volume();
            const double Vj = shapes[j].volume();
            const double ri = Dimension::rootnu(Vi/M_PI);
            const double rj = Dimension::rootnu(Vj/M_PI);
            const Vector dji = centers[j] - centers[i];
            const double djimag = dji.magnitude();
            const Vector delta = dji.unitVector() * dispfrac*max(0.1*min(ri, rj), max(0.0, ri + rj - djimag));
            if (flags[j] == 3) {
              displacements[i] -= 0.25*Vi/(Vi + Vj)*delta;
              displacements[j] += 0.25*Vj/(Vi + Vj)*delta;
            } else {
              displacements[i] -= 0.25*delta;
            }
          }
        }
      }
    }

    // Apply the displacements of this iteration.
    for (auto i = imin; i < imax; ++i) {
      if (flags[i] >= 2) {
        centers[i] += displacements[i];
      }
    }

#ifdef USE_MPI
    // Global broadcast of the new centers.
    for (int iproc = 0; iproc != nprocs; ++iproc) {
      int jmin = imin, jmax = imax;
      MPI_Bcast(&jmin, 1, MPI_INT, iproc, Communicator::communicator());
      MPI_Bcast(&jmax, 1, MPI_INT, iproc, Communicator::communicator());
      vector<char> buffer;
      int bufsize = 0;
      if (rank == iproc) {
        for (auto i = imin; i != imax; ++i) packElement(centers[i], buffer);
        bufsize = buffer.size();
      }
      MPI_Bcast(&bufsize, 1, MPI_INT, iproc, Communicator::communicator());
      buffer.resize(bufsize);
      MPI_Bcast(&buffer[0], bufsize, MPI_CHAR, iproc, Communicator::communicator());
      if (rank != iproc) {
        vector<char>::const_iterator bufitr = buffer.begin();
        for (auto j = jmin; j < jmax; ++j) unpackElement(centers[j], bufitr, buffer.end());
        CHECK(bufitr == buffer.end());
      }
    }
#endif

    //  Check the current level of overlap.
    maxoverlap = 0.0;
    for (auto i = imin; i < imax; ++i) {
      if (flags[i] >= 2) {
        flags[i] = 2;
        const FacetedVolume shapei = shapes[i] + centers[i];
        const double Ri = radius[shape2node[i]];
        const vector<vector<int>>& neighbors = cm.connectivityForNode(0, shape2node[i]);
        CHECK(neighbors.size() == 1);
        for (auto j: neighbors[0]) {
          CHECK(j != shape2node[i]);
          if (flags[j] >= 1) {
            FacetedVolume shapej;
            if (flags[j] == 1) {
              shapej = shapes[j];
            } else {
              shapej = shapes[j] + centers[j];
            }
            if (shapei.intersect(shapej)) {
                const Vector centj = shapej.centroid();
                const double Rj = radius[shape2node[j]];
                auto overlap = max(0.0, (Ri + Rj - (shapei.centroid() - shapej.centroid()).magnitude())/(Ri + Rj));

                // # Ri = (shapei.volume/pi)**(1.0/self.ndim)
                // # Rj = (shapej.volume/pi)**(1.0/self.ndim)
                // # vol = min(shapei.volume, shapej.volume)
                // # verticesi = shapei.vertices()
                // # verticesj = shapej.vertices()
                // # dmax = 0.0
                // # for vi in verticesi:
                // #     if shapej.contains(vi):
                // #         dmax = max(dmax, shapej.distance(vi))
                // # for vj in verticesj:
                // #     if shapei.contains(vj):
                // #         dmax = max(dmax, shapei.distance(vj))
                // # overlap = dmax/(vol/pi)**(1.0/self.ndim)

                if (overlap > maxoverlapfrac) {
                  flags[i] = 3;
                  maxoverlap = max(maxoverlap, overlap);
                }
              }
          }
        }
      }
    }

#ifdef USE_MPI
    // Global broadcast of the new flags.
    for (int iproc = 0; iproc != nprocs; ++iproc) {
      int jmin = imin, jmax = imax;
      MPI_Bcast(&jmin, 1, MPI_INT, iproc, Communicator::communicator());
      MPI_Bcast(&jmax, 1, MPI_INT, iproc, Communicator::communicator());
      vector<char> buffer;
      int bufsize = 0;
      if (rank == iproc) {
        for (auto i = imin; i != imax; ++i) packElement(flags[i], buffer);
        bufsize = buffer.size();
      }
      MPI_Bcast(&bufsize, 1, MPI_INT, iproc, Communicator::communicator());
      buffer.resize(bufsize);
      MPI_Bcast(&buffer[0], bufsize, MPI_CHAR, iproc, Communicator::communicator());
      if (rank != iproc) {
        vector<char>::const_iterator bufitr = buffer.begin();
        for (auto j = jmin; j < jmax; ++j) unpackElement(flags[j], bufitr, buffer.end());
        CHECK(bufitr == buffer.end());
      }
    }
    maxoverlap = allReduce(maxoverlap, MPI_MAX, Communicator::communicator());
    // print "   Iteration %i, maxoverlap %g" % (iter, maxoverlap)
#endif

    // Any shapes we were unable to disentangle turn back to inactive, otherwise set the successful
    // survivors to flag=1.
    for (auto i = imin; i < imax; ++i) {
      if (flags[i] == 2) {
        flags[i] = 1;
      } else if (flags[i] == 3) {
        CHECK(maxoverlap > maxoverlapfrac);
        flags[i] = 0;
        // # # Make one last ditch attempt to randomly fit this shape in.
        // # centers[i] = self.randomCenter(i, centers)
        // # if centers[i]:
        // #     flags[i] = 1
        // # else:
        // #     flags[i] = 0
      }
    }

#ifdef USE_MPI
    // Global broadcast of the new geometry.
    for (int iproc = 0; iproc != nprocs; ++iproc) {
      int jmin = imin, jmax = imax;
      MPI_Bcast(&jmin, 1, MPI_INT, iproc, Communicator::communicator());
      MPI_Bcast(&jmax, 1, MPI_INT, iproc, Communicator::communicator());
      vector<char> buffer;
      int bufsize = 0;
      if (rank == iproc) {
        for (auto i = imin; i != imax; ++i) packElement(flags[i], buffer);
        for (auto i = imin; i != imax; ++i) packElement(centers[i], buffer);
        bufsize = buffer.size();
      }
      MPI_Bcast(&bufsize, 1, MPI_INT, iproc, Communicator::communicator());
      buffer.resize(bufsize);
      MPI_Bcast(&buffer[0], bufsize, MPI_CHAR, iproc, Communicator::communicator());
      if (rank != iproc) {
        vector<char>::const_iterator bufitr = buffer.begin();
        for (auto j = jmin; j < jmax; ++j) unpackElement(flags[j], bufitr, buffer.end());
        for (auto j = jmin; j < jmax; ++j) unpackElement(centers[j], bufitr, buffer.end());
        CHECK(bufitr == buffer.end());
      }
    }
#endif

  } // end of iteration

  // That's it.
  return iter;
}

}
