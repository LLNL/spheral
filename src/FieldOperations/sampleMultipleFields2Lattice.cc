//---------------------------------Spheral++----------------------------------//
// sampleMultipleFields2Lattice.
//
// SPH sample all the Fields in a FieldListSet to a lattice of positions.
// The results are returned as a 
// tuple< vector< vector<Scalar> >,
//        vector< vector<Vector> >,
//        vector< vector<Tensor> >,
//        vector< vector<SymTensor> >
//
// Created by JMO, Wed Nov 16 10:40:07 PST 2005
//----------------------------------------------------------------------------//
#include <algorithm>

#include "sampleMultipleFields2Lattice.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

#include "DBC.hh"

#ifdef USE_MPI
#include "mpi.h"
#include "Distributed/BoundingVolumeDistributedBoundary.hh"
#include "Distributed/Communicator.hh"
#endif

namespace Spheral {
namespace FieldSpace {

using namespace std;

using boost::tuple;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using KernelSpace::TableKernel;

//------------------------------------------------------------------------------
// Compute the step size.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
stepSize(const typename Dimension::Vector& xmin,
         const typename Dimension::Vector& xmax,
         const vector<int>& nsample) {
  REQUIRE(nsample.size() == Dimension::nDim);
  for (int i = 0; i != Dimension::nDim; ++i) REQUIRE(nsample[i] > 0);

  typename Dimension::Vector result = xmax - xmin;
  for (int i = 0; i != Dimension::nDim; ++i) result(i) /= nsample[i];
  return result;
}

//------------------------------------------------------------------------------
// Generate the dimension appropriate index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
vector<int>
makeIndex(const int ix,
          const int iy,
          const int iz) {
  REQUIRE(false);
}

template<>
inline
vector<int>
makeIndex<Dim<1> >(const int ix,
                   const int iy,
                   const int iz) {
  vector<int> result(1);
  result[0] = ix;
  return result;
}

template<>
inline
vector<int>
makeIndex<Dim<2> >(const int ix,
                   const int iy,
                   const int iz) {
  vector<int> result(2);
  result[0] = ix;
  result[1] = iy;
  return result;
}

template<>
inline
vector<int>
makeIndex<Dim<3> >(const int ix,
                   const int iy,
                   const int iz) {
  vector<int> result(3);
  result[0] = ix;
  result[1] = iy;
  result[2] = iz;
  return result;
}

//------------------------------------------------------------------------------
// Compute the indicies of the lattice overlapping the given point.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
vector< vector<int> >
latticePoints(const typename Dimension::Vector& ri,
              const typename Dimension::Vector& extent,
              const typename Dimension::Vector& xmin,
              const typename Dimension::Vector& stepSize,
              const vector<int>& nsample) {

  // Pre-conditions.
  REQUIRE(extent.minElement() > 0.0);
  REQUIRE(stepSize.minElement() > 0.0);
  REQUIRE(nsample.size() == Dimension::nDim);
  for (int i = 0; i != Dimension::nDim; ++i) REQUIRE(nsample[i] > 0);

  // Find the min & max indicies in each dimension.
  vector<int> imin(size_t(3), 0);
  vector<int> imax(size_t(3), 0);
  int ntot = 1;
  for (int i = 0; i != Dimension::nDim; ++i) {
    imin[i] = max(0, min(nsample[i] - 1, int((ri(i) - extent(i) - xmin(i))/stepSize(i))));
    imax[i] = max(0, min(nsample[i] - 1, int((ri(i) + extent(i) - xmin(i))/stepSize(i)) + 1));
    CHECK(imax[i] - imin[i] >= 0);
    ntot *= imax[i] - imin[i] + 1;
  }
  CHECK(ntot >= 0);

  // Now build the result with the full set of indicies.
  vector< vector<int> > result;
  result.reserve(ntot);
  int ix = imin[0];
  int iy = imin[1];
  int iz = imin[2];
  int i = 0;
  while (i < ntot) {
    CHECK(ix >= imin[0] && ix <= imax[0]);
    CHECK(iy >= imin[1] && iy <= imax[1]);
    CHECK(iz >= imin[2] && iz <= imax[2]);
    result.push_back(makeIndex<Dimension>(ix, iy, iz));
    ++i;
    ++ix;
    if (ix > imax[0]) {
      ix = imin[0];
      ++iy;
      if (iy > imax[1]) {
        iy = imin[1];
        ++iz;
      }
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  ENSURE(result.size() == ntot);
  for (int i = 0; i != result.size(); ++i) {
    ENSURE(result[i].size() == Dimension::nDim);
    for (int j = 0; j != Dimension::nDim; ++j) {
      ENSURE(result[i][j] >= imin[j] && result[i][j] <= imax[j]);
    }
  }
  END_CONTRACT_SCOPE;

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Compute the position for the given lattice index.
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
latticePosition(const vector<int>& indicies,
                const Dim<1>::Vector& xmin,
                const Dim<1>::Vector& xmax,
                const Dim<1>::Vector& xstep) {
  REQUIRE(indicies.size() == 1);
  const Dim<1>::Vector result(xmin.x() + indicies[0]*xstep.x());
  REQUIRE(result >= xmin && result <= xmax);
  return result;
}

inline
Dim<2>::Vector
latticePosition(const vector<int>& indicies,
                const Dim<2>::Vector& xmin,
                const Dim<2>::Vector& xmax,
                const Dim<2>::Vector& xstep) {
  REQUIRE(indicies.size() == 2);
  const Dim<2>::Vector result(xmin.x() + indicies[0]*xstep.x(),
                              xmin.y() + indicies[1]*xstep.y());
  REQUIRE(result >= xmin && result <= xmax);
  return result;
}

inline
Dim<3>::Vector
latticePosition(const vector<int>& indicies,
                const Dim<3>::Vector& xmin,
                const Dim<3>::Vector& xmax,
                const Dim<3>::Vector& xstep) {
  REQUIRE(indicies.size() == 3);
  const Dim<3>::Vector result(xmin.x() + indicies[0]*xstep.x(),
                              xmin.y() + indicies[1]*xstep.y(),
                              xmin.z() + indicies[2]*xstep.z());
  REQUIRE(result >= xmin && result <= xmax);
  return result;
}

//------------------------------------------------------------------------------
// A simple minded mapping of lattice index to domain.
//------------------------------------------------------------------------------
inline
int
domainForLatticeIndex(const int i,
                      const int nlocal,
                      const int remainder) {
  REQUIRE(nlocal > 0);
  const int off = remainder*(nlocal + 1);
  return min(i, off)/(nlocal + 1) + max(0, i - off)/nlocal;
}

//------------------------------------------------------------------------------
// Convert the tuple of indicies into the flat index for the lattice arrays.
//------------------------------------------------------------------------------
inline
int
latticeIndex(const vector<int>& indicies,
             const vector<int>& nsample) {
  REQUIRE(indicies.size() == nsample.size());
  REQUIRE(nsample.size() == 1 ||
          nsample.size() == 2 ||
          nsample.size() == 3);
  if (nsample.size() == 1) {
    return indicies[0];
  } else if (nsample.size() == 2) {
    return indicies[1]*nsample[0] + indicies[0];
  } else {
    return (indicies[2]*nsample[1]*nsample[0] + 
            indicies[1]*nsample[0] + 
            indicies[0]);
  }
}

inline
int
localLatticeIndex(const int iglobal,
                  const int nlocal,
                  const int remainder) {
  const int domain = domainForLatticeIndex(iglobal, nlocal, remainder);
  return iglobal - min(domain, remainder)*(nlocal + 1) - max(0, domain - remainder)*nlocal;
}

//------------------------------------------------------------------------------
// The sample function itself.
//------------------------------------------------------------------------------
template<typename Dimension>
void
sampleMultipleFields2Lattice(const FieldListSet<Dimension>& fieldListSet,
                             const FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldList<Dimension, typename Dimension::Scalar>& weight,
                             const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                             const FieldList<Dimension, int>& mask,
                             const TableKernel<Dimension>& W,
                             const typename Dimension::Vector& xmin,
                             const typename Dimension::Vector& xmax,
                             const vector<int>& nsample,
                             vector< vector<typename Dimension::Scalar> >& scalarValues,
                             vector< vector<typename Dimension::Vector> >& vectorValues,
                             vector< vector<typename Dimension::Tensor> >& tensorValues,
                             vector< vector<typename Dimension::SymTensor> >& symTensorValues) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Sizes of stuff.
  const size_t numScalarFieldLists = fieldListSet.ScalarFieldLists.size();
  const size_t numVectorFieldLists = fieldListSet.VectorFieldLists.size();
  const size_t numTensorFieldLists = fieldListSet.TensorFieldLists.size();
  const size_t numSymTensorFieldLists = fieldListSet.SymTensorFieldLists.size();

  // Pre-conditions.
  for (typename vector< FieldList<Dimension, Scalar> >::const_iterator fieldListItr = fieldListSet.ScalarFieldLists.begin();
       fieldListItr != fieldListSet.ScalarFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Scalar>& fieldList = *fieldListItr;
    for (int i = 0; i != fieldList.numFields(); ++i) {
      VERIFY(position.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(weight.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(Hfield.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(mask.haveNodeList(*fieldList[i]->nodeListPtr()));
    }
  }
  for (typename vector< FieldList<Dimension, Vector> >::const_iterator fieldListItr = fieldListSet.VectorFieldLists.begin();
       fieldListItr != fieldListSet.VectorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Vector>& fieldList = *fieldListItr;
    for (int i = 0; i < fieldList.numFields(); ++i) {
      VERIFY(position.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(weight.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(Hfield.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(mask.haveNodeList(*fieldList[i]->nodeListPtr()));
    }
  }
  for (typename vector< FieldList<Dimension, Tensor> >::const_iterator fieldListItr = fieldListSet.TensorFieldLists.begin();
       fieldListItr != fieldListSet.TensorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Tensor>& fieldList = *fieldListItr;
    for (int i = 0; i < fieldList.numFields(); ++i) {
      VERIFY(position.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(weight.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(Hfield.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(mask.haveNodeList(*fieldList[i]->nodeListPtr()));
    }
  }
  for (typename vector< FieldList<Dimension, SymTensor> >::const_iterator fieldListItr = fieldListSet.SymTensorFieldLists.begin();
       fieldListItr != fieldListSet.SymTensorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, SymTensor>& fieldList = *fieldListItr;
    for (int i = 0; i < fieldList.numFields(); ++i) {
      VERIFY(position.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(weight.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(Hfield.haveNodeList(*fieldList[i]->nodeListPtr()));
      VERIFY(mask.haveNodeList(*fieldList[i]->nodeListPtr()));
    }
  }
  VERIFY(position.numFields() == weight.numFields());
  VERIFY(position.numFields() == Hfield.numFields());
  VERIFY(position.numFields() == mask.numFields());
  for (int i = 0; i != position.numFields(); ++i) {
    VERIFY(position[i]->nodeListPtr() == weight[i]->nodeListPtr());
    VERIFY(position[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    VERIFY(position[i]->nodeListPtr() == mask[i]->nodeListPtr());
  }
  VERIFY(xmin < xmax);
  VERIFY(nsample.size() == Dimension::nDim);
  for (int i = 0; i != Dimension::nDim; ++i) VERIFY(nsample[i] > 0);

  // Get the parallel geometry.
  int procID = 0;
  int numProcs = 1;
#ifdef USE_MPI
  MPI_Comm_rank(Communicator::communicator(), &procID);
  MPI_Comm_size(Communicator::communicator(), &numProcs);
#endif
  CHECK(numProcs > 0);

  // We need to exclude any nodes that come from the Distributed boundary condition.
#ifdef USE_MPI
  BoundarySpace::BoundingVolumeDistributedBoundary<Dimension>& distributedBoundary = BoundarySpace::BoundingVolumeDistributedBoundary<Dimension>::instance();
#endif

  // Compute the total number of sample points.
  int ntotal = 1;
  for (int i = 0; i != Dimension::nDim; ++i) ntotal *= nsample[i];
  CHECK(ntotal > 0);

  // Compute the local number of sample points (needed as input to domainForLatticeIndex.
  const int nlocal = ntotal / numProcs;
  const int nremainder = ntotal % numProcs;
  CHECK(nlocal*numProcs + nremainder == ntotal);

  // Compute the step size.
  const Vector xstep = stepSize<Dimension>(xmin, xmax, nsample);

  // Size of the kernel.
  const double kernelExtent = W.kernelExtent();

  // Prepare the result.
  const int nlocalsizing = nlocal + (procID < nremainder ? 1 : 0);
  scalarValues = vector<vector<Scalar> >();
  vectorValues = vector<vector<Vector> >();
  tensorValues = vector<vector<Tensor> >();
  symTensorValues = vector<vector<SymTensor> >();
  for (int i = 0; i != numScalarFieldLists; ++i) scalarValues.push_back(vector<Scalar>(nlocalsizing, 0.0));
  for (int i = 0; i != numVectorFieldLists; ++i) vectorValues.push_back(vector<Vector>(nlocalsizing, Vector()));
  for (int i = 0; i != numTensorFieldLists; ++i) tensorValues.push_back(vector<Tensor>(nlocalsizing, Tensor()));
  for (int i = 0; i != numSymTensorFieldLists; ++i) symTensorValues.push_back(vector<SymTensor>(nlocalsizing, SymTensor()));

  // This data structure is how we store the locally computed values, which
  // are then reduced in parallel to the final result.
  // Note this is indexed by the global lattice index, not local.
  typedef map<int, tuple<vector<Scalar>,
                         vector<Vector>,
                         vector<Tensor>,
                         vector<SymTensor> > > LocalStorage;
  typedef tuple<vector<Scalar>,
                vector<Vector>,
                vector<Tensor>,
                vector<SymTensor> > LocalElement;
  LocalStorage localResult;

  // Loop over the positions.
  for (AllNodeIterator<Dimension> nodeItr = position.nodeBegin();
       nodeItr != position.nodeEnd();
       ++nodeItr) {

#ifdef USE_MPI
    const bool useNode = count(distributedBoundary.ghostNodes(*(nodeItr.nodeListPtr())).begin(),
                               distributedBoundary.ghostNodes(*(nodeItr.nodeListPtr())).end(),
                               nodeItr.nodeID()) == 0;
#else
    const bool useNode = true;
#endif
    if (useNode and mask(nodeItr) == 1) {

      // Sample node (i) state.
      const Vector& ri = position(nodeItr);
      const SymTensor& Hi = Hfield(nodeItr);
      const Scalar& weighti = weight(nodeItr);
      const Vector extent = nodeItr.nodeListPtr()->neighbor().HExtent(Hi, kernelExtent);
      const Scalar Hdeti = Hi.Determinant();

      // Compute the set of lattice positions this node contributes to.
      const vector< vector<int> > ids = latticePoints<Dimension>(ri, extent, 
                                                                 xmin, xstep, nsample);

      // Iterate over the points.
      for (vector< vector<int> >::const_iterator itr = ids.begin();
           itr != ids.end();
           ++itr) {

        // Compute the sample position and index into the lattice sample array.
        const Vector rj = latticePosition(*itr, xmin, xmax, xstep);
        const int j = latticeIndex(*itr, nsample);
        CHECK(rj >= xmin && rj <= xmax);
        CHECK(j >= 0 && j < ntotal);

        // Do we have an entry for this point?
        if (localResult.find(j) == localResult.end()) {
          localResult[j] = LocalElement(vector<Scalar>(numScalarFieldLists),
                                        vector<Vector>(numVectorFieldLists),
                                        vector<Tensor>(numTensorFieldLists),
                                        vector<SymTensor>(numSymTensorFieldLists));
        }

        // Contribution of the SPH point to this position.
        const Vector rij = ri - rj;
        const Scalar etai = (Hi*rij).magnitude();
        const Scalar Wi = W(etai, Hdeti);
        const Scalar thpt = weighti*Wi;

        // Scalar fields.
        {
          vector<Scalar>& samples = boost::tuples::get<0>(localResult[j]);
          for (int k = 0; k != numScalarFieldLists; ++k) {
            CHECK(k < samples.size());
            const FieldList<Dimension, Scalar>& fieldList = fieldListSet.ScalarFieldLists[k];
            if (fieldList.haveNodeList(*nodeItr.nodeListPtr()))
              samples[k] += thpt*fieldList(nodeItr);
          }
        }
          
        // Vector fields.
        {
          vector<Vector>& samples = boost::tuples::get<1>(localResult[j]);
          for (int k = 0; k != numVectorFieldLists; ++k) {
            CHECK(k < samples.size());
            const FieldList<Dimension, Vector>& fieldList = fieldListSet.VectorFieldLists[k];
            if (fieldList.haveNodeList(*nodeItr.nodeListPtr()))
              samples[k] += thpt*fieldList(nodeItr);
          }
        }
          
        // Tensor fields.
        {
          vector<Tensor>& samples = boost::tuples::get<2>(localResult[j]);
          for (int k = 0; k != numTensorFieldLists; ++k) {
            CHECK(k < samples.size());
            const FieldList<Dimension, Tensor>& fieldList = fieldListSet.TensorFieldLists[k];
            if (fieldList.haveNodeList(*nodeItr.nodeListPtr()))
              samples[k] += thpt*fieldList(nodeItr);
          }
        }
          
        // SymTensor fields.
        {
          vector<SymTensor>& samples = boost::tuples::get<3>(localResult[j]);
          for (int k = 0; k != numSymTensorFieldLists; ++k) {
            CHECK(k < samples.size());
            const FieldList<Dimension, SymTensor>& fieldList = fieldListSet.SymTensorFieldLists[k];
            if (fieldList.haveNodeList(*nodeItr.nodeListPtr()))
              samples[k] += thpt*fieldList(nodeItr);
          }
        }
          
      }
    }
  }

  // In parallel we have to reduce the elements across processors.

  // Calculate the size of the packed data per position.
  const int sizeOfElement = (numScalarFieldLists*DataTypeTraits<Scalar>::numElements(0.0)*sizeof(typename DataTypeTraits<Scalar>::ElementType) +
                             numVectorFieldLists*DataTypeTraits<Vector>::numElements(Vector::zero)*sizeof(typename DataTypeTraits<Vector>::ElementType) +
                             numTensorFieldLists*DataTypeTraits<Tensor>::numElements(Tensor::zero)*sizeof(typename DataTypeTraits<Tensor>::ElementType) +
                             numSymTensorFieldLists*DataTypeTraits<SymTensor>::numElements(SymTensor::zero)*sizeof(typename DataTypeTraits<SymTensor>::ElementType));

  // Figure out what we have to send to other processors.
  // In the process we transfer any of our local values we've accumulated to the final result.
  vector< vector<int> > sendIndiciesBuffers(numProcs);
  vector< vector<char> > sendValuesBuffers(numProcs);
  for (typename LocalStorage::const_iterator itr = localResult.begin();
       itr != localResult.end();
       ++itr) {
    const int j = itr->first;
    const int jdomain = domainForLatticeIndex(j, nlocal, nremainder);
    const int jlocal = localLatticeIndex(j, nlocal, nremainder);

    // Does this element belong to this processor?
    if (jdomain == procID) {
      CHECK(jlocal >= 0 and jlocal < nlocalsizing);
      
      // Scalar fields.
      {
        const vector<Scalar>& localSamples = boost::tuples::get<0>(itr->second);
        for (int k = 0; k != numScalarFieldLists; ++k) {
          CHECK(k < localSamples.size());
          CHECK(k < scalarValues.size() and jlocal < scalarValues[k].size());
          scalarValues[k][jlocal] += localSamples[k];
        }
      }
      
      // Vector fields.
      {
        const vector<Vector>& localSamples = boost::tuples::get<1>(itr->second);
        for (int k = 0; k != numVectorFieldLists; ++k) {
          CHECK(k < localSamples.size());
          CHECK(k < vectorValues.size() and jlocal < vectorValues[k].size());
          vectorValues[k][jlocal] += localSamples[k];
        }
      }
      
      // Tensor fields.
      {
        const vector<Tensor>& localSamples = boost::tuples::get<2>(itr->second);
        for (int k = 0; k != numTensorFieldLists; ++k) {
          CHECK(k < localSamples.size());
          CHECK(k < tensorValues.size() and jlocal < tensorValues[k].size());
          tensorValues[k][jlocal] += localSamples[k];
        }
      }
      
      // SymTensor fields.
      {
        const vector<SymTensor>& localSamples = boost::tuples::get<3>(itr->second);
        for (int k = 0; k != numSymTensorFieldLists; ++k) {
          CHECK(k < localSamples.size());
          CHECK(k < symTensorValues.size() and jlocal < symTensorValues[k].size());
          symTensorValues[k][jlocal] += localSamples[k];
        }
      }

    } else {

      // We have to send this data to another processor, so pack it up.
      CHECK(jdomain < sendIndiciesBuffers.size());
      CHECK(jdomain < sendValuesBuffers.size());
      sendIndiciesBuffers[jdomain].push_back(jlocal);

      // Scalar fields.
      {
        const vector<Scalar>& localSamples = boost::tuples::get<0>(itr->second);
        for (int k = 0; k != numScalarFieldLists; ++k) {
          CHECK(k < localSamples.size());
          packElement(localSamples[k], sendValuesBuffers[jdomain]);
        }
      }

      // Vector fields.
      {
        const vector<Vector>& localSamples = boost::tuples::get<1>(itr->second);
        for (int k = 0; k != numVectorFieldLists; ++k) {
          CHECK(k < localSamples.size());
          packElement(localSamples[k], sendValuesBuffers[jdomain]);
        }
      }

      // Tensor fields.
      {
        const vector<Tensor>& localSamples = boost::tuples::get<2>(itr->second);
        for (int k = 0; k != numTensorFieldLists; ++k) {
          CHECK(k < localSamples.size());
          packElement(localSamples[k], sendValuesBuffers[jdomain]);
        }
      }

      // SymTensor fields.
      {
        const vector<SymTensor>& localSamples = boost::tuples::get<3>(itr->second);
        for (int k = 0; k != numSymTensorFieldLists; ++k) {
          CHECK(k < localSamples.size());
          packElement(localSamples[k], sendValuesBuffers[jdomain]);
        }
      }

    }
  }

#ifdef USE_MPI
  // Send everyone all the information we have for them.
  vector<int> numSends(numProcs);
  vector<MPI_Request> sendRequests(3*(numProcs - 1));
  for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
    if (sendProc != procID) {
      const int bufIndex = sendProc > procID ? sendProc - 1 : sendProc;
      numSends[sendProc] = sendIndiciesBuffers[sendProc].size();
      CHECK(sendValuesBuffers[sendProc].size() == numSends[sendProc]*sizeOfElement);
      MPI_Isend(&numSends[sendProc], 1, MPI_INT, sendProc, 1, Communicator::communicator(), &sendRequests[bufIndex]);
      MPI_Isend(&(*sendIndiciesBuffers[sendProc].begin()), numSends[sendProc], MPI_INT, sendProc, 2, Communicator::communicator(), &sendRequests[(numProcs - 1) + bufIndex]);
      MPI_Isend(&(*sendValuesBuffers[sendProc].begin()), numSends[sendProc]*sizeOfElement, MPI_CHAR, sendProc, 3, Communicator::communicator(), &sendRequests[2*(numProcs - 1) + bufIndex]);
    }
  }

  // Post receives to see how many indicies other processors are sending to us.
  vector<int> numReceiveNodes(size_t(numProcs), 0);
  vector<MPI_Request> recvRequests0(numProcs - 1);
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    if (recvProc != procID) {
      const int bufIndex = recvProc > procID ? recvProc - 1 : recvProc;
      MPI_Irecv(&numReceiveNodes[recvProc], 1, MPI_INT, recvProc, 1, Communicator::communicator(), &recvRequests0[bufIndex]);
    }
  }

  // Wait until we have the sizes from everyone.
  {
    vector<MPI_Status> recvStatus(recvRequests0.size());
    MPI_Waitall(recvRequests0.size(), &(*recvRequests0.begin()), &(*recvStatus.begin()));
  }

  // Size up the receive buffers.
  vector<vector<int>  > recvIndiciesBuffers(numProcs);
  vector<vector<char> > recvValuesBuffers(numProcs);
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    recvIndiciesBuffers[recvProc].resize(numReceiveNodes[recvProc]);
    recvValuesBuffers[recvProc].resize(numReceiveNodes[recvProc]*sizeOfElement);
  }

  // Post receives for the nodal data.
  vector<MPI_Request> recvRequests1(2*(numProcs - 1));
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    if (recvProc != procID) {
      const int bufIndex = recvProc > procID ? recvProc - 1 : recvProc;
      MPI_Irecv(&(*recvIndiciesBuffers[recvProc].begin()), numReceiveNodes[recvProc], MPI_INT, recvProc, 2, Communicator::communicator(), &recvRequests1[bufIndex]);
      MPI_Irecv(&(*recvValuesBuffers[recvProc].begin()), numReceiveNodes[recvProc]*sizeOfElement, MPI_CHAR, recvProc, 3, Communicator::communicator(), &recvRequests1[numProcs - 1 + bufIndex]);
    }
  }

  // Wait until we have the full receive data.
  {
    vector<MPI_Status> recvStatus(recvRequests1.size());
    MPI_Waitall(recvRequests1.size(), &(*recvRequests1.begin()), &(*recvStatus.begin()));
  }

  // Now unpack the receive values and add them onto our local values.
  for (int recvProc = 0; recvProc != numProcs; ++recvProc) {
    if (recvProc != procID) {
      const vector<char>& buffer = recvValuesBuffers[recvProc];
      CHECK(buffer.size() % sizeOfElement == 0);
      vector<char>::const_iterator bufItr = buffer.begin();
      for (int i = 0; i != numReceiveNodes[recvProc]; ++i) {
        const int jlocal = recvIndiciesBuffers[recvProc][i];
        CHECK(jlocal >= 0 and jlocal < nlocalsizing);

        // Scalar Fields.
        {
          Scalar element;
          for (int k = 0; k != numScalarFieldLists; ++k) {
            CHECK(k < scalarValues.size() and jlocal < scalarValues[k].size());
            unpackElement(element, bufItr, buffer.end());
            scalarValues[k][jlocal] += element;
          }
        }

        // Vector Fields.
        {
          Vector element;
          for (int k = 0; k != numVectorFieldLists; ++k) {
            CHECK(k < vectorValues.size() and jlocal < vectorValues[k].size());
            unpackElement(element, bufItr, buffer.end());
            vectorValues[k][jlocal] += element;
          }
        }

        // Tensor Fields.
        {
          Tensor element;
          for (int k = 0; k != numTensorFieldLists; ++k) {
            CHECK(k < tensorValues.size() and jlocal < tensorValues[k].size());
            unpackElement(element, bufItr, buffer.end());
            tensorValues[k][jlocal] += element;
          }
        }

        // SymTensor Fields.
        {
          SymTensor element;
          for (int k = 0; k != numSymTensorFieldLists; ++k) {
            CHECK(k < symTensorValues.size() and jlocal < symTensorValues[k].size());
            unpackElement(element, bufItr, buffer.end());
            symTensorValues[k][jlocal] += element;
          }
        }

      }
      CHECK(bufItr == buffer.end());
    }
  }

  // Wait until all our sends are completed.
  {
    vector<MPI_Status> sendStatus(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &(*sendRequests.begin()), &(*sendStatus.begin()));
  }

#endif

  // That's it.
}

}
}

