//----------------------------------------------------------------------------
//! \author $Author: jeffjohnson $
//! \date $Date: 2002-11-05 23:24:26 -0800 (Tue, 05 Nov 2002) $
//! \version $Revision: 461 $
//
// Implementation of the LaplacianMatrixFactory class.
//----------------------------------------------------------------------------

namespace Spheral {

//----------------------------------------------------------------------------
template <typename Dimension>
double** 
LaplacianMatrixFactory<Dimension>::
SPHGradSquared(const DataBase<Dimension>& db)
{
  // Count up some particles.
  size_t numberOfParticles = FIXME;
  
  // First, we need to allocate storage for the grad and gradSquared matrices.
  // They should be N x N, where N is the number of particles.
  double** grad = new double*[numberOfParticles];
  double** gradSquared = new double*[numberOfParticles];
  for (size_t i = 0; i < numberOfParticles)
  {
    grad[i] = new double[numberOfParticles];
    gradSquared[i] = new double[numberOfParticles];
    std::fill(gradSquared[i], gradSquared[i] + numberOfParticles, 0.0);
  } // end for

  // Now we've got to build the damn grad operator.

  // Some important field lists.
  FieldList<Dimension, Vector>& position = db.globalPosition();
  
  // Loop over each particle...
  for (typename DataBase<Dimension>::IDIterator
       ithNodeIter = db.internalNodeBegin();
       ithNodeIter = db.internalNodeEnd();
       ++ithNodeIter)
  {
    int i = ithNodeIter.nodeID();
    const Vector& ri = position(ithNodeIter);
    const SymTensor& Hi = Hfield(ithNodeIter);

    // Set the refined neighbor information for this master node.
    dataBase.setRefineNeighborNodeLists(position(masterItr), Hfield(masterItr));
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator nodeListItr = 
        dataBase.fluidNodeListBegin();
        nodeListItr < dataBase.fluidNodeListEnd();
        ++nodeListItr) {
      updateCoarseNeighborStats((*nodeListItr)->neighborPtr()->numCoarse());
      updateRefineNeighborStats((*nodeListItr)->neighborPtr()->numRefine());
    } // end for

    // Loop over the neighbors of this particle.
    for (IDIterator jthNodeIter = db.refineNodeBegin();
         jthNodeIter < db.refineNodeEnd();
         ++jthNodeIter)
    {
      // This neighbor will be particle j.
      j = jthNodeIter.nodeID();
      rj = position(jthNodeIter);
      Hj = Hfield(jthNodeIter);

      // Calculate rij, the displacement between particles i and j.
      Vector rij = ri - rj;

      // Calculate the normalized displacements.
      Vector etai = Hi * rij;
      Vector etaiNorm = etai.unitVector();
      Vector etaj = Hj * rij;
      Vector etajNorm = etaj.unitVector();

      grad[i][j] = kernel.grad(r, h);
    } // end for

  } // end for

  // Now we square the grad matrix to obtain gradSquared.
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    for (size_t j = 0; j < numberOfParticles; ++j)
    {
      for (size_t k = 0; k < numberOfParticles; ++k)
      {
        gradSquared[i][j] += grad[i][k] * grad[k][j];
      } // end for
    } // end for
  } // end for

  // Toss the grad matrix.
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    delete [] grad[i];
  } // end for
  delete [] grad;
  
  // Now return the gradSquared matrix.
  return gradSquared;
} // end SPHGradSquared
//----------------------------------------------------------------------------

} // end Spheral

