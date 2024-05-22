//---------------------------------Spheral++----------------------------------//
// Resistive SPMHD implementation.
//
//! \author $Author: jeffjohnson $
//! \version $Revision: 2239 $
//! \date $Date: 2007-05-28 23:58:39 -0700 (Mon, 28 May 2007) $
//----------------------------------------------------------------------------//
#include <algorithm>

#include "MHD/MHD.hh"
#include "MHD/MHDFieldNames.hh"
#include "MHD/ConductingFluidNodeList.hh"
#include "MHD/CurrentDensityUpdatePolicy.hh"
#include "MHD/MagnetosonicSpeedUpdatePolicy.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/NonDynamicState.hh"
#include "Hydro/WeightPolicy.hh"
#include "Hydro/OmegaGradhPolicy.hh"
#include "Material/EquationOfState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Boundary/PeriodicBoundary.hh"
#include "Distributed/DistributedBoundary.hh"
#include "Utilities/DBC.hh"
#include "Material/PhysicalConstants.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Spasmos/Config.h"
#include "Spasmos/PyPetsc.h"
#include "Spasmos/MatFactory.h"
#include "Utilities/Communicator.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>

namespace Spheral {


//------------------------------------------------------------------------------
MHD::
MHD(const TableKernel<Dim<3> >& kernel,
    double mu0,
    double implicitness):
  mDiffResNorm(0.0),
  mDiffNumIters(0),
  mMu0(mu0),
  mImplicitness(0.5),
  mKernel(kernel),
  mVecFactory(0),
  mMatFactory(0),
  mDiffMatrix(0),
  mDiffRHS(0),
  mIsPeriodicOrDistributedNode(),
  mDiffSolver(),
  mNodeIndices(FieldList<Dim<3>, int>::Copy),
  mOverlapNodes(),
  mMagneticEnergy(0.),
  mStabilizer(maxStressStabilizer), 
  mS0(),
  mBext(),
  mDivBCleaner(hyperbolicCleaner),
  mSpecificThermalEnergyUpdate(integrateSpecificThermalEnergy),
  mTotalSpecificEnergyUpdate(dontUpdateTotalSpecificEnergy),
  mPsi(FieldList<Dim<3>, Dim<3>::Scalar>::Copy),
  mDpsiDt(FieldList<Dim<3>, Dim<3>::Scalar>::Copy),
  mMaxDivB(-FLT_MAX),
  mMinDivB(FLT_MAX),
  mAvgDivB(0.0),
  mFirstStep(true) { 


  // Make sure Spasmos is around.
  importConfig();
  Spasmos_Initialize(SPASMOS_ERRORHANDLER_EXCEPTION);
  importPyPetsc();

  // Initialize matrix and vector factories.
  mVecFactory = VecFactory_New();
  mMatFactory = MatFactory_New();

  // Initialize a GMRES solver for the magnetic diffusion equation.
  KSPCreate(PETSC_COMM_WORLD, &mDiffSolver);
  KSPSetType(mDiffSolver, KSPGMRES);

  // Use a Jacobi preconditioner.
  PC precond;
  KSPGetPC(mDiffSolver, &precond);
  PCSetType(precond, PCJACOBI);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
MHD::
~MHD() {
   // Get rid of the diffusion matrix.
   Py_XDECREF(mDiffMatrix);

   // Kill the linear solver.
   KSPDestroy(mDiffSolver);

   // Kill the vector and matrix factories.
   Py_DECREF(mVecFactory);
   Py_DECREF(mMatFactory);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void 
MHD::
mComputeMatrixStructure(const DataBase<Dim<3> >& dataBase,
                        const State<Dim<3> >& state) const
{
  // Initialize the matrix with a nonzero structure.

  // Refresh our global indexing scheme for all the nodes in the database.
  mNodeIndices = globalNodeIDs(dataBase);
  int n = mNodeIndices.numInternalNodes();

  // Make sure that the ghost nodes have the right global indexing.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mNodeIndices);
  }

  // This is necessary to ease the pain of the parallel boundary condition.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->finalizeGhostBoundary();
  }

  // Additionally, we mark ghost nodes associated with parallel and periodic 
  // boundary conditions as being special.
  mIsPeriodicOrDistributedNode.clear();
  const vector<FluidNodeList<Dim<3> >*>& nodeLists = dataBase.fluidNodeListPtrs();
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    if ((dynamic_cast<DistributedBoundary<Dim<3> >*>(*boundaryItr) != 0) or 
        (dynamic_cast<PeriodicBoundary<Dim<3> >*>(*boundaryItr) != 0))
    {
       // Get the ghost nodes for this Boundary for each node list.
       for (size_t i = 0; i < nodeLists.size(); ++i)
       {
          const vector<int>& ghostNodes = (*boundaryItr)->ghostNodes(*(nodeLists[i]));
          for (size_t j = 0; j < ghostNodes.size(); ++j)
          {
             // Mark this node as special.
             pair<int, int> nodeEntry(i, ghostNodes[j]);
             mIsPeriodicOrDistributedNode[nodeEntry] = 1;
          }
       }
    }
  }

  // Now create an "adjacency table" that contains the nonzero columns in each 
  // row for the diffusion matrix.
  vector<vector<Py_ssize_t> > adjTable(3*n);
  const ConnectivityMap<Dim<3> >& connectivityMap = dataBase.connectivityMap();
  int numNodeLists = static_cast<int>(nodeLists.size());
  for (size_t iNodeList = 0; iNodeList < numNodeLists; ++iNodeList)
  {
    // Traverse the internal nodes in this node list.
    for (size_t iNode = 0; iNode < nodeLists[iNodeList]->numInternalNodes(); ++iNode)
    {
      // Keep track of the rows and columns in the non-zero block we create for 
      // the ith node, and automatically include the diagonal (ii) element.
      Py_ssize_t iNodeIndex = (*mNodeIndices[iNodeList])[iNode];
      Py_ssize_t row = iNodeIndex;
      adjTable[3*iNode].push_back(3*iNodeIndex);
      adjTable[3*iNode+1].push_back(3*iNodeIndex+1);
      adjTable[3*iNode+2].push_back(3*iNodeIndex+2);
       
      // Get the neighbors for this node.
      const vector<vector<int> >& jNodes = 
        connectivityMap.connectivityForNode(nodeLists[iNodeList], iNode);

      // Iterate over the nodelists containing the ith node's neighbors.
      for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
      {
        // Get the neighboring nodes within this node list.
        const vector<int>& neighbors = jNodes[jNodeList];

        // Get the number of internal nodes within the jth nodelist.
        int numInternalNodes = nodeLists[jNodeList]->numInternalNodes();

        // Traverse these neighboring nodes.
        for (vector<int>::const_iterator jNodeIter = neighbors.begin();
             jNodeIter != neighbors.end(); ++jNodeIter)
        {
          int jNode = *jNodeIter;

          // Scrutinize non-internal nodes.  If a non-internal node is a 
          // boundary node related to a periodic or a distributed boundary, 
          // we add its contribution to the matrix.  Otherwise we discard it.
          if (jNode >= numInternalNodes)
          {
            pair<int, int> nodeEntry(jNodeList, jNode);
            if (!mIsPeriodicOrDistributedNode[nodeEntry])
              continue;
          }

          // Add the jth node to the list of rows and columns in our block.
          int jNodeIndex = ((*(mNodeIndices[jNodeList]))[jNode]);
#if 0
if ((jNodeIndex - (*mNodeIndices[0])[0] >= mNodeIndices.numInternalNodes()) ||
    (jNodeIndex < (*mNodeIndices[0])[0]))
   printf("%d: %d is off-proc\n", iNodeIndex, jNodeIndex);
#endif
          adjTable[3*iNode].push_back(static_cast<Py_ssize_t>(3*jNodeIndex));
          adjTable[3*iNode+1].push_back(static_cast<Py_ssize_t>(3*jNodeIndex+1));
          adjTable[3*iNode+2].push_back(static_cast<Py_ssize_t>(3*jNodeIndex+2));
        }
      }
    }
  }

  // Build the diffusion matrix from the given non-zero structure.
  Py_XDECREF(mDiffMatrix);
  mDiffMatrix = MatFactory_SparseMatrix_FromTable(mMatFactory, adjTable);

  // Tell the matrix not to allow insertions outside of the prescribed non-zeros.
  // FIXME: Something is wrong!  I don't think the nonzeros are set properly 
  // FIXME: in the parallel case!
//  MatSetOption(PyPetsc_MAT(mDiffMatrix), MAT_NEW_NONZERO_LOCATION_ERR);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
Mat 
MHD::
mCurlCurlMatrix(const DataBase<Dim<3> >& dataBase,
                const State<Dim<3> >& state) const
{
  REQUIRE(mDiffMatrix != 0);

  // Get the FieldLists storing the data we will use in computing our matrix.
  const FieldList<Dim<3>, Scalar> weight = state.scalarFields(HydroFieldNames::weight);
  const FieldList<Dim<3>, Scalar> resistivity = state.scalarFields(MHDFieldNames::resistivity);
  const FieldList<Dim<3>, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dim<3>, SymTensor> Htensor = state.symTensorFields(HydroFieldNames::H);
  const vector<FluidNodeList<Dim<3> >*>& nodeLists = dataBase.fluidNodeListPtrs();
  int numNodeLists = static_cast<int>(nodeLists.size());
  CHECK(weight.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(Htensor.size() == numNodeLists);

  // Create a matrix in which to store the curl-curl operator.  Since mDiffMatrix
  // exists and holds the non-zero structure that we desire, we simply duplicate 
  // it without copying the non-zero values.
  Mat D;
  MatDuplicate(PyPetsc_MAT(mDiffMatrix), MAT_DO_NOT_COPY_VALUES, &D);

  // Clear any overlapping node entries.
  mOverlapNodes.clear();

  // Traverse the list of nodes in our database and add their contributions 
  // to the matrix.
  const ConnectivityMap<Dim<3> >& connectivityMap = dataBase.connectivityMap();
  for (size_t iNodeList = 0; iNodeList < numNodeLists; ++iNodeList)
  {
    const Field<Dim<3>, Scalar>& weighti = *weight[iNodeList];
    const Field<Dim<3>, Scalar>& Ri = *resistivity[iNodeList];
    const Field<Dim<3>, Vector>& positioni = *position[iNodeList];
    const Field<Dim<3>, SymTensor>& Htensori = *Htensor[iNodeList];

    // Traverse the internal nodes in this node list.
    for (size_t iNode = 0; iNode < nodeLists[iNodeList]->numInternalNodes(); ++iNode)
    {
      // If the ith node is one that overlaps another node, we skip its 
      // row in the matrix.
      int nodeIndex = (*(mNodeIndices[iNodeList]))[iNode];
      PetscInt rows[3];
      rows[0] = 3*nodeIndex;
      rows[1] = 3*nodeIndex+1;
      rows[2] = 3*nodeIndex+2;
      if (mOverlapNodes.find(nodeIndex) != mOverlapNodes.end())
         continue;

      // Get data for the ith node.
      Scalar wi = weighti[iNode];
      Scalar Di = Ri[iNode]/mMu0;
      Vector xi = positioni[iNode];
      SymTensor Hi = Htensori[iNode];
      Scalar detHi = Hi.Determinant();

      // Get the neighbors for this node and count them up.
      size_t numNeighbors = 0;
      const vector<vector<int> >& jNodes = 
        connectivityMap.connectivityForNode(nodeLists[iNodeList], iNode);
      for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
         numNeighbors += jNodes[jNodeList].size();

      // Keep track of the rows and columns in the non-zero block we create for 
      // the ith node, and include the diagonal (ii) element.
      vector<PetscInt> columns(3*(1+numNeighbors));
      vector<double> values(3*(1+numNeighbors));
      
      // The diagonal entry should be the first one in the row, regardless of how 
      // things are actually ordered.
      columns[0] = rows[0];
      columns[1] = rows[1];
      columns[2] = rows[2];
      values[0] = values[1] = values[2] = 0.0;
       
      // If our resistivity isn't spatially homogeneous, we need to compute its
      // gradient so that we can include its effects on our matrix.
      Vector gradDi;
      State<Dim<3> >::ScalarPolicyPointerType RPolicy = state.policyForField(Ri);
      if (dynamic_cast<NonDynamicState<Dim<3>, Field<Dim<3>, Scalar> >*>(RPolicy.get()) == 0)
      {
        // Our policy is not a NonDynamicState, so we have a resistivity that 
        // may be spatially inhomogeneous.  Let's compute its gradient.
        for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
        {
          // Get the neighboring nodes within this node list.
          const vector<int>& neighbors = jNodes[jNodeList];
          const Field<Dim<3>, Scalar>& weightj = *weight[jNodeList];
          const Field<Dim<3>, Scalar>& Rj = *resistivity[jNodeList];
          const Field<Dim<3>, Vector>& positionj = *position[jNodeList];
          const Field<Dim<3>, SymTensor>& Htensorj = *Htensor[jNodeList];

          // Traverse these neighboring nodes.
          for (vector<int>::const_iterator jNodeIter = neighbors.begin();
               jNodeIter != neighbors.end(); ++jNodeIter)
          {
            // Get data for the node.
            int jNode = *jNodeIter;
            Vector xj = positionj[jNode];
            Scalar wj = weightj[jNode];
            Scalar Dj = Rj[jNode]/mMu0;
            SymTensor Hj = Htensorj[jNode];
            Scalar detHj = Hj.Determinant();

            // Node displacement.
            const Vector xij = xi - xj;
            const Vector etai = Hi*xij;
            const Vector etaj = Hj*xij;
            const Scalar etaMagi = etai.magnitude();
            const Scalar etaMagj = etaj.magnitude();
            CHECK(etaMagi >= 0.0);
            CHECK(etaMagj >= 0.0);

            // Symmetrized kernel weight and gradient.
            const Scalar Wi = mKernel(etaMagi, detHi);
            const Vector Hetai = Hi*etai.unitVector();
            const Vector gradWi = Hetai*mKernel.grad(etaMagi, detHi);

            const Scalar Wj = mKernel(etaMagj, detHj);
            const Vector Hetaj = Hj*etaj.unitVector();
            const Vector gradWj = Hetaj*mKernel.grad(etaMagj, detHj);

            const Vector gradWij = 0.5 * (gradWi + gradWj);

            // Gradient contribution.
            gradDi += wj * (Dj - Di) * gradWij;
          }
        }
      }

      // Iterate over the nodelists containing the ith node's neighbors and 
      // determine the matrix elements associated with the contributing 
      // neighbors.
      for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
      {
        // Get the neighboring nodes within this node list.
        const vector<int>& neighbors = jNodes[jNodeList];
        const Field<Dim<3>, Scalar>& weightj = *weight[jNodeList];
        const Field<Dim<3>, Scalar>& Rj = *resistivity[jNodeList];
        const Field<Dim<3>, Vector>& positionj = *position[jNodeList];
        const Field<Dim<3>, SymTensor>& Htensorj = *Htensor[jNodeList];

        // Traverse these neighboring nodes.
        for (size_t j = 0; j < neighbors.size(); ++j)
        {
          // Scrutinize non-internal nodes.  If a non-internal node is a 
          // boundary node related to a periodic or a distributed boundary, 
          // we add its contribution to the matrix.  Otherwise we discard it.
          int jNode = neighbors[j];
          if (jNode >= nodeLists[jNodeList]->numInternalNodes())
          {
            pair<int, int> nodeEntry(jNodeList, jNode);
            if (!mIsPeriodicOrDistributedNode[nodeEntry])
              continue;
          }

          // Get data for the node.
          Vector xj = positionj[jNode];
          Scalar wj = weightj[jNode];
          Scalar Dj = Rj[jNode]/mMu0;
          SymTensor Hj = Htensorj[jNode];
          Scalar detHj = Hj.Determinant();

          // Find the global index for this neighbor node.
          PetscInt jNodeIndex = 
            static_cast<PetscInt>((*(mNodeIndices[jNodeList]))[jNode]);

          // Node displacement.
          const Vector xij = xi - xj;
          const Vector etai = Hi*xij;
          const Vector etaj = Hj*xij;
          const Scalar etaMagi = etai.magnitude();
          const Scalar etaMagj = etaj.magnitude();
          CHECK(etaMagi >= 0.0);
          CHECK(etaMagj >= 0.0);

          // Does node j overlap node i?  That is, are the two so close 
          // together that they are physically indistinguishable?  If so,
          // skip this entry of the matrix.
          if (etaMagi < 0.1)
          {
             mOverlapNodes[jNodeIndex] = nodeIndex;
//printf("Overlapping node %d -> %d\n", jNodeIndex, row);
             continue;
          }

          // Add the jth node to the list of columns.
          columns.push_back(jNodeIndex);

          // Symmetrized kernel weight and gradient.
          const Scalar Wi = mKernel(etaMagi, detHi);
          const Vector Hetai = Hi*etai.unitVector();
          const Vector gradWi = Hetai*mKernel.grad(etaMagi, detHi);

          const Scalar Wj = mKernel(etaMagj, detHj);
          const Vector Hetaj = Hj*etaj.unitVector();
          const Vector gradWj = Hetaj*mKernel.grad(etaMagj, detHj);

          const Vector gradWij = 0.5 * (gradWi + gradWj);

          // Now compute the off-diagonal value, softening the interaction 
          // of nodes that are close together.
          double hi = pow(detHi, 1.0/3.0);
          double Rij2 = xij.magnitude2() + 0.01*hi*hi;
          double offDiagContrib[3];
          offDiagContrib[0] = offDiagContrib[1] = offDiagContrib[2] = 
            2.0 * Di * wj * gradWij.dot(xij)/Rij2 // constant resistivity terms 
            - wj * gradDi.dot(gradWij);           // resistivity gradient terms
          offDiagContrib[0] += wj * gradDi.x()*gradWij.x();
          offDiagContrib[1] += wj * gradDi.y()*gradWij.y();
          offDiagContrib[2] += wj * gradDi.z()*gradWij.z();
          values[3*j] += offDiagContrib[0];
          values[3*numNeighbors + 3*j + 1] += offDiagContrib[1];
          values[6*numNeighbors + 3*j + 2] += offDiagContrib[2];

          // The diagonal value for this row is the negative sum of the off-diagonal
          // values, so we simply subtract.
          values[0] -= offDiagContrib[0];
          values[1] -= offDiagContrib[1];
          values[2] -= offDiagContrib[2];
        }
      }
      assert(values.size() == columns.size());

      // Now insert the row into the matrix. 
      if (MatSetValues(D, 3, &(rows[0]), 
                       static_cast<PetscInt>(columns.size()), &(columns[0]),
                       &(values[0]), INSERT_VALUES) != 0)
      {
         PyErr_Print();
         VERIFY(false);
      }
    }
  }

  // Before we leave, we must assemble the matrix.
  MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);

  // Return our matrix.
  return D;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
MHD::
mApplyDiffusionBCs(const DataBase<Dim<3> >& dataBase, 
                   const State<Dim<3> >& state,
                   PyObject* matrix, 
                   PyObject* RHS) const
{
  // Go over the set of overlapping nodes and mash together the corresponding 
  // values in the right-hand side vector.
  if (!mOverlapNodes.empty())
  {
    vector<PetscInt> overlappingRows, degenerateRows;
    for (map<int, int>::const_iterator iter = mOverlapNodes.begin();
         iter != mOverlapNodes.end(); ++iter)
    {
      degenerateRows.push_back(static_cast<PetscInt>(iter->first));
      overlappingRows.push_back(static_cast<PetscInt>(iter->second));
    }
    vector<double> degenerateDensities(mOverlapNodes.size());
    VecGetValues(PyPetsc_VEC(RHS), static_cast<PetscInt>(degenerateRows.size()),
          &(degenerateRows[0]), &(degenerateDensities[0]));
    VecSetValues(PyPetsc_VEC(RHS), static_cast<PetscInt>(overlappingRows.size()),
          &(overlappingRows[0]), &(degenerateDensities[0]), ADD_VALUES);
    VecAssemblyBegin(PyPetsc_VEC(RHS));
    VecAssemblyEnd(PyPetsc_VEC(RHS));

    // Make sure to zero out the entries on the degenerate rows, too. 
    fill(degenerateDensities.begin(), degenerateDensities.end(), 0.0);
    VecSetValues(PyPetsc_VEC(RHS), static_cast<PetscInt>(degenerateRows.size()),
          &(degenerateRows[0]), &(degenerateDensities[0]), INSERT_VALUES);
    VecAssemblyBegin(PyPetsc_VEC(RHS));
    VecAssemblyEnd(PyPetsc_VEC(RHS));
  }

  // Add contributions from boundary nodes to the RHS.
  const FieldList<Dim<3>, Scalar> weight = state.scalarFields(HydroFieldNames::weight);
  const FieldList<Dim<3>, Scalar> resistivity = state.scalarFields(MHDFieldNames::resistivity);
  const FieldList<Dim<3>, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dim<3>, Vector> induction = state.vectorFields(MHDFieldNames::magneticInduction);
  const FieldList<Dim<3>, SymTensor> Htensor = state.symTensorFields(HydroFieldNames::H);
  const vector<FluidNodeList<Dim<3> >*>& nodeLists = dataBase.fluidNodeListPtrs();
  const ConnectivityMap<Dim<3> >& connectivityMap = dataBase.connectivityMap();
  size_t numNodeLists = nodeLists.size();
  for (size_t iNodeList = 0; iNodeList < dataBase.numFluidNodeLists(); ++iNodeList)
  {
    const Field<Dim<3>, Scalar>& weighti = *weight[iNodeList];
    const Field<Dim<3>, Scalar>& Ri = *resistivity[iNodeList];
    const Field<Dim<3>, Vector>& positioni = *position[iNodeList];
    const Field<Dim<3>, SymTensor>& Htensori = *Htensor[iNodeList];

    // Traverse the internal nodes in this node list.
    for (size_t iNode = 0; iNode < nodeLists[iNodeList]->numInternalNodes(); ++iNode)
    {
      // If the ith node is one that overlaps another node, we skip its 
      // row in the matrix.
      int nodeIndex = (*(mNodeIndices[iNodeList]))[iNode];
      PetscInt rows[3];
      rows[0] = 3*nodeIndex;
      rows[1] = 3*nodeIndex+1;
      rows[2] = 3*nodeIndex+2;
      if (mOverlapNodes.find(nodeIndex) != mOverlapNodes.end())
         continue;

      // Get data for the ith node.
      Scalar wi = weighti[iNode];
      Scalar Di = Ri[iNode]/mMu0;
      Vector xi = positioni[iNode];
      SymTensor Hi = Htensori[iNode];
      Scalar detHi = Hi.Determinant();

      // Get the neighbors for this node and count them up.
      size_t numNeighbors = 0;
      const vector<vector<int> >& jNodes = 
        connectivityMap.connectivityForNode(nodeLists[iNodeList], iNode);
      for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
         numNeighbors += jNodes[jNodeList].size();

      // Here are the values we'll insert into the vector.
      double values[3] = {0.0, 0.0, 0.0};
       
      // If our resistivity isn't spatially homogeneous, we need to compute its
      // gradient so that we can include its effects on our matrix.
      Vector gradDi;
      State<Dim<3> >::ScalarPolicyPointerType RPolicy = state.policyForField(Ri);
      if (dynamic_cast<NonDynamicState<Dim<3>, Field<Dim<3>, Scalar> >*>(RPolicy.get()) == 0)
      {
        // Our policy is not a NonDynamicState, so we have a resistivity that 
        // may be spatially inhomogeneous.  Let's compute its gradient.
        for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
        {
          // Get the neighboring nodes within this node list.
          const vector<int>& neighbors = jNodes[jNodeList];
          const Field<Dim<3>, Scalar>& weightj = *weight[jNodeList];
          const Field<Dim<3>, Scalar>& Rj = *resistivity[jNodeList];
          const Field<Dim<3>, Vector>& positionj = *position[jNodeList];
          const Field<Dim<3>, SymTensor>& Htensorj = *Htensor[jNodeList];

          // Traverse these neighboring nodes.
          for (vector<int>::const_iterator jNodeIter = neighbors.begin();
               jNodeIter != neighbors.end(); ++jNodeIter)
          {
            // Get data for the node.
            int jNode = *jNodeIter;
            Vector xj = positionj[jNode];
            Scalar wj = weightj[jNode];
            Scalar Dj = Rj[jNode]/mMu0;
            SymTensor Hj = Htensorj[jNode];
            Scalar detHj = Hj.Determinant();

            // Node displacement.
            const Vector xij = xi - xj;
            const Vector etai = Hi*xij;
            const Vector etaj = Hj*xij;
            const Scalar etaMagi = etai.magnitude();
            const Scalar etaMagj = etaj.magnitude();
            CHECK(etaMagi >= 0.0);
            CHECK(etaMagj >= 0.0);

            // Symmetrized kernel weight and gradient.
            const Scalar Wi = mKernel(etaMagi, detHi);
            const Vector Hetai = Hi*etai.unitVector();
            const Vector gradWi = Hetai*mKernel.grad(etaMagi, detHi);

            const Scalar Wj = mKernel(etaMagj, detHj);
            const Vector Hetaj = Hj*etaj.unitVector();
            const Vector gradWj = Hetaj*mKernel.grad(etaMagj, detHj);

            const Vector gradWij = 0.5 * (gradWi + gradWj);

            // Gradient contribution.
            gradDi += wj * (Dj - Di) * gradWij;
          }
        }
      }

      // Iterate over the nodelists containing the ith node's neighbors and 
      // determine the matrix elements associated with the contributing 
      // neighbors.
      for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
      {
        // Get the neighboring nodes within this node list.
        const vector<int>& neighbors = jNodes[jNodeList];
        const Field<Dim<3>, Scalar>& weightj = *weight[jNodeList];
        const Field<Dim<3>, Scalar>& Rj = *resistivity[jNodeList];
        const Field<Dim<3>, Vector>& positionj = *position[jNodeList];
        const Field<Dim<3>, Vector>& inductionj = *induction[jNodeList];
        const Field<Dim<3>, SymTensor>& Htensorj = *Htensor[jNodeList];

        // Traverse these neighboring nodes.
        for (size_t j = 0; j < neighbors.size(); ++j)
        {
          // Skip internal nodes and periodic/distributed boundary nodes. 
          int jNode = neighbors[j];
          if (jNode >= nodeLists[jNodeList]->numInternalNodes())
          {
            pair<int, int> nodeEntry(jNodeList, jNode);
            if (mIsPeriodicOrDistributedNode[nodeEntry])
              continue;
          }
          else
             continue; // Skip internal nodes.

          // Get data for the node.
          Vector xj = positionj[jNode];
          Vector Bj = inductionj[jNode];
          Scalar wj = weightj[jNode];
          Scalar Dj = Rj[jNode]/mMu0;
          SymTensor Hj = Htensorj[jNode];
          Scalar detHj = Hj.Determinant();

          // Find the global index for this neighbor node.
          PetscInt jNodeIndex = 
            static_cast<PetscInt>((*(mNodeIndices[jNodeList]))[jNode]);

          // Node displacement.
          const Vector xij = xi - xj;
          const Vector etai = Hi*xij;
          const Vector etaj = Hj*xij;
          const Scalar etaMagi = etai.magnitude();
          const Scalar etaMagj = etaj.magnitude();
          CHECK(etaMagi >= 0.0);
          CHECK(etaMagj >= 0.0);

          // Symmetrized kernel weight and gradient.
          const Scalar Wi = mKernel(etaMagi, detHi);
          const Vector Hetai = Hi*etai.unitVector();
          const Vector gradWi = Hetai*mKernel.grad(etaMagi, detHi);

          const Scalar Wj = mKernel(etaMagj, detHj);
          const Vector Hetaj = Hj*etaj.unitVector();
          const Vector gradWj = Hetaj*mKernel.grad(etaMagj, detHj);

          const Vector gradWij = 0.5 * (gradWi + gradWj);

          // Now compute the off-diagonal value, softening the interaction 
          // of nodes that are close together.
          double hi = pow(detHi, 1.0/3.0);
          double Rij2 = xij.magnitude2() + 0.01*hi*hi;
          Vector jContrib = -2.0*Di*wj*Bj*gradWij.dot(xij)/Rij2 + 
                            wj*Bj*gradDi.dot(gradWij) - 
                            wj*Bj.dyad(gradWij).dot(gradDi);
          values[0] += jContrib.x();
          values[1] += jContrib.y();
          values[2] += jContrib.z();
        }
      }

      // Now insert the values into the vector.
      if (VecSetValues(PyPetsc_VEC(RHS), 3, &(rows[0]), &(values[0]), INSERT_VALUES) != 0)
      {
         PyErr_Print();
         VERIFY(false);
      }
    }
  }

  // Assemble the vector.
  if ((VecAssemblyBegin(PyPetsc_VEC(RHS)) != 0) ||
      (VecAssemblyEnd(PyPetsc_VEC(RHS)) != 0))
    VERIFY(false);

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
MHD::
mAddDiffusionDerivatives(const DataBase<Dim<3> >& dataBase,
                         const State<Dim<3> >& state,
                         const Dim<3>::Scalar time,
                         const Dim<3>::Scalar dt,
                         StateDerivatives<Dim<3> >& derivs) const {
  // First, make sure that the matrix is up to date. 
  REQUIRE(mDiffMatrix != 0);

  // First we compute the curl-curl matrix representing the magnetic diffusion 
  // term in the induction evolution equation.
  Mat D = mCurlCurlMatrix(dataBase, state);

  // Now form the left-hand side of the linear system (which we call the diffusion
  // matrix): LHS = 1 + alpha*dt*D.
  MatCopy(D, PyPetsc_MAT(mDiffMatrix), SAME_NONZERO_PATTERN);
  MatScale(PyPetsc_MAT(mDiffMatrix), mImplicitness*dt);
  MatShift(PyPetsc_MAT(mDiffMatrix), 1.0);

  // Create a vector representing B at time tn.
  Py_ssize_t n = static_cast<Py_ssize_t>(mNodeIndices.numInternalNodes());
  PyObject* Bn = VecFactory_Vector(mVecFactory, 3*n);
  CHECK(Bn != 0);
  const FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);
  vector<double> localVec(3*B.numInternalNodes());
  size_t offset = 0;
  for (int iNodeList = 0; iNodeList < B.numFields(); ++iNodeList) {
    Field<Dim<3>, Vector>& Bi = *(B[iNodeList]);
    int numInternalNodes = Bi.nodeList().numInternalNodes();
    for (int iNode = 0; iNode < numInternalNodes; ++iNode) {
      localVec[3*offset] = Bi[iNode].x();
      localVec[3*offset+1] = Bi[iNode].y();
      localVec[3*offset+2] = Bi[iNode].z();
      ++offset;
    }
  }
  CHECK(offset == B.numInternalNodes());
  PyVec_CopyFromArray(Bn, 0, static_cast<Py_ssize_t>(localVec.size()),
                      &(localVec[0]));

  // Form the right-hand side of the linear system by multiplying Bn by
  // (1-implicitness)*dt*D.
  Py_XDECREF(mDiffRHS);
  mDiffRHS = VecFactory_Vector(mVecFactory, 3*n);
  CHECK(mDiffRHS != 0);
  MatMult(D, PyPetsc_VEC(Bn), PyPetsc_VEC(mDiffRHS));
  VecAXPBY(PyPetsc_VEC(mDiffRHS), 1.0, -dt*(1.0 - mImplicitness), 
           PyPetsc_VEC(Bn));

  // At this point, we're done with D.
  MatDestroy(D);

  // Apply any pertinent boundary conditions.
  mApplyDiffusionBCs(dataBase, state, mDiffMatrix, mDiffRHS);

  // Now set up the linear system and solve it.
  PyObject* Bn1 = VecFactory_Vector(mVecFactory, 3*n);
  KSPSetOperators(mDiffSolver, PyPetsc_MAT(mDiffMatrix), 
                  PyPetsc_MAT(mDiffMatrix), SAME_PRECONDITIONER);
  KSPSolve(mDiffSolver, PyPetsc_VEC(mDiffRHS), PyPetsc_VEC(Bn1));

  // Get data about the linear solution.
  KSPGetResidualNorm(mDiffSolver, &mDiffResNorm);
  KSPGetIterationNumber(mDiffSolver, &mDiffNumIters);

  // Did the solution converge?  If not, we have to let the user know What 
  // Happened.  The easiest way to do this is to let Spasmos decipher the 
  // PETSc error.
  if (!PyPetsc_SolveConverged((void*)mDiffSolver, 
                              "SPH magnetic diffusion equation"))
  {
    PyErr_Print();
    VERIFY(false);
  }

  // Now go over the set of overlapping nodes and copy the solutions back 
  // to the associated degenerate nodes.
  // together.
  if (!mOverlapNodes.empty())
  {
    vector<PetscInt> overlappingRows, degenerateRows;
    for (map<int, int>::const_iterator iter = mOverlapNodes.begin();
         iter != mOverlapNodes.end(); ++iter)
    {
      degenerateRows.push_back(static_cast<PetscInt>(iter->first));
      overlappingRows.push_back(static_cast<PetscInt>(iter->second));
//printf("Copying solution from %d to %d\n", iter->second, iter->first);
    }
    vector<double> degenerateDensities(mOverlapNodes.size());
    VecGetValues(PyPetsc_VEC(Bn1), static_cast<PetscInt>(overlappingRows.size()),
                 &(overlappingRows[0]), &(degenerateDensities[0]));
    VecSetValues(PyPetsc_VEC(Bn1), static_cast<PetscInt>(degenerateRows.size()),
                 &(degenerateRows[0]), &(degenerateDensities[0]), INSERT_VALUES);
    VecAssemblyBegin(PyPetsc_VEC(Bn1));
    VecAssemblyEnd(PyPetsc_VEC(Bn1));
  }

  // Compute dB/dt by dividing the change in the magnetic induction by the 
  // timestep, storing the result in Bn1. 
  VecAXPBY(PyPetsc_VEC(Bn1), -1.0/dt, 1.0/dt, PyPetsc_VEC(Bn));

  // Now we're done with Bn.
  Py_DECREF(Bn);

  // Add the diffusive dB/dt terms to the time derivatives for B.
  FieldList<Dim<3>, Vector> DBDt = derivs.vectorFields(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::magneticInduction);
  PyVec_CopyToArray(Bn1, 0, static_cast<Py_ssize_t>(localVec.size()),
                    &(localVec[0]));
  for (int iNodeList = 0; iNodeList < mNodeIndices.numFields(); ++iNodeList) {
    Field<Dim<3>, Vector>& DBiDt = *(DBDt[iNodeList]);
    int numInternalNodes = DBiDt.nodeList().numInternalNodes();
    for (int iNode = 0; iNode < numInternalNodes; ++iNode) {
      DBiDt[iNode] += Vector(localVec[3*iNode], localVec[3*iNode+1], localVec[3*iNode+2]);
    }
  }

  // Now we're done with Bn1!
  Py_DECREF(Bn1);

  // And in fact we're done, period.
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// This helper function determines whether a set of materials has resistivity.
bool 
hasResistivity(const DataBase<Dim<3> >& dataBase,
               const State<Dim<3> >& state)
{
  // We don't assume there's resistivity till we see some evidence of it.
  bool haveResistivity = false;

  FieldList<Dim<3>, Dim<3>::Scalar> R = state.scalarFields(MHDFieldNames::resistivity);
  const ConnectivityMap<Dim<3> >& connectivityMap = dataBase.connectivityMap();
  int numNodeLists = R.numFields();
  for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
    const NodeList<Dim<3> >& nodeList = R[iNodeList]->nodeList(); 

    // If the nodelist is empty, we define it to be a perfect conductor
    // (for NOW).
    if (nodeList.numNodes() == 0) continue;

    // Is this a conducting fluid?  If not, we can skip it, since there's no 
    // magnetic field otherwise.
    const ConductingFluidNodeList* condFluidPtr = 
      dynamic_cast<const ConductingFluidNodeList*>(&nodeList);
    if (condFluidPtr == 0) continue;

    // Check the resistivity update policy to see whether there's any 
    // resistivity.  If the resistivity model is non-dynamic and the value of 
    // the resistivity is 0, then it's a perfect conductor.  Otherwise, 
    // there is resistivity.
    if (dynamic_cast<const NonDynamicState<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >*>(condFluidPtr->resistivityPolicy().get()) != 0) {
      const Field<Dim<3>, Dim<3>::Scalar>& Ri = condFluidPtr->resistivity();
      double minR = *min_element(Ri.begin(), Ri.end());
      double maxR = *max_element(Ri.begin(), Ri.end());
      if ((minR != 0.0) || (maxR != 0.0))
      {
        haveResistivity = true;
      }
    }
    else {
      haveResistivity = true;
    }
  }
  return haveResistivity;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void 
MHD::
evaluateDerivatives(const Dim<3>::Scalar time,
                    const Dim<3>::Scalar dt,
                    const DataBase<Dim<3> >& dataBase,
                    const State<Dim<3> >& state,
                    StateDerivatives<Dim<3> >& derivs) const
{

  // Access to pertinent fields in the database.
  const FieldList<Dim<3>, Scalar> m = state.scalarFields(HydroFieldNames::mass);
  const FieldList<Dim<3>, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dim<3>, Vector> x = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dim<3>, Vector> v = state.vectorFields(HydroFieldNames::velocity);
  const FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);
  const FieldList<Dim<3>, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  const FieldList<Dim<3>, Scalar> Omega = state.scalarFields(HydroFieldNames::omegaGradh);
  const FieldList<Dim<3>, Scalar> u = state.scalarFields(HydroFieldNames::specificThermalEnergy);
  const FieldList<Dim<3>, Scalar> P = state.scalarFields(HydroFieldNames::pressure);
  const FieldList<Dim<3>, Scalar> psi = state.scalarFields(MHDFieldNames::hyperbolicCleaning);
  const FieldList<Dim<3>, Scalar> divB = state.scalarFields(MHDFieldNames::magneticDivergence);

  // Get the accelerations and inductions we'll be modifying.
  FieldList<Dim<3>, Vector> DvDt = derivs.vectorFields(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + HydroFieldNames::velocity);
  FieldList<Dim<3>, Vector> DBDt = derivs.vectorFields(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::magneticInduction);
  FieldList<Dim<3>, Scalar> DeDt = derivs.scalarFields(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::totalSpecificEnergy);
  FieldList<Dim<3>, Scalar> DpsiDt = derivs.scalarFields(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::hyperbolicCleaning);
  FieldList<Dim<3>, vector<Vector> > pairAccelerations = derivs.vectorVectorFields(HydroFieldNames::pairAccelerations);

  // Should we use resistive MHD?
  bool isResistive = hasResistivity(dataBase, state);

  // Now compute the time derivatives of v, B, and possibly e.
  const ConnectivityMap<Dim<3> >& connectivityMap = dataBase.connectivityMap();
  int numNodeLists = DvDt.numFields();
double maxAcc = 0.0;
  for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
    const NodeList<Dim<3> >& nodeList = DvDt[iNodeList]->nodeList(); 

    // Is this a conducting fluid?  If not, we can skip it, since there's no 
    // magnetic field otherwise.
    const ConductingFluidNodeList* condFluidPtr = 
      dynamic_cast<const ConductingFluidNodeList*>(&nodeList);
    if (condFluidPtr == 0) continue;

    // Now compute those derivatives.
    for (int iNode = 0; iNode < nodeList.numInternalNodes(); ++iNode) {
      // Get a reference to the acceleration vector and the induction.
      Scalar mi = (*m[iNodeList])[iNode];
      Vector& dVdti = (*DvDt[iNodeList])[iNode];
      const Vector& Bi = (*B[iNodeList])[iNode];
      Vector& dBdti = (*DBDt[iNodeList])[iNode];
      const SymTensor& Hi = (*H[iNodeList])[iNode];
      const Vector& Vi = (*v[iNodeList])[iNode];
      Scalar detHi = Hi.Determinant();
      Scalar rhoi = (*rho[iNodeList])[iNode];
      double Pi = (*P[iNodeList])[iNode];
      double Omegai = (*Omega[iNodeList])[iNode];
      vector<Vector>& pairAcci = (*pairAccelerations[iNodeList])[iNode];

      // Potential (and gradient) used for hyperbolic divergence cleaning.
      double psii = 0.0;
      if (psi.numFields() > 0)
         psii = (*psi[iNodeList])[iNode];
      Vector gradPsi;

      // Compute the Maxwell stress tensor.
      Scalar Bi2 = Bi.magnitude2();
      SymTensor Si = (1.0/mMu0) * (Bi.selfdyad() - 
                                   0.5 * Bi2 * SymTensor(1., 0., 0.,
                                                         0., 1., 0.,
                                                         0., 0., 1.)) - mS0;
      
      // Loop over the neighboring nodes.
      const vector< vector<int> >& fullConnectivity = 
         connectivityMap.connectivityForNode(&nodeList, iNode);
      // Iterate over the NodeLists.
      for (int jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
        const vector<int>& connectivity = fullConnectivity[jNodeList];
        if (connectivity.size() > 0) {
          // Loop over the neighbors.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            int jNode = *jItr;
            CHECK(iNode != jNode);
            Scalar mj = (*m[jNodeList])[jNode];
            Scalar rhoj = (*rho[jNodeList])[jNode];
            double Pj = (*P[jNodeList])[jNode];
            const Vector& Bj = (*B[jNodeList])[jNode];
            Vector xij = (*x[iNodeList])[iNode] - (*x[jNodeList])[jNode];
            CHECK(xij.magnitude2() != 0.0);
            const Vector& Vj = (*v[jNodeList])[jNode];
            Vector Vij = Vi - Vj;
            double Omegaj = (*Omega[jNodeList])[jNode];
            double psij = 0.0;
            if (psi.numFields() > 0)
               psij = (*psi[jNodeList])[jNode];
            vector<Vector>& pairAccj = (*pairAccelerations[jNodeList])[jNode];

            const SymTensor& Hj = (*H[jNodeList])[jNode];
            Scalar detHj = Hj.Determinant();

            // Compute the symmetrized kernel.
            Vector etai = Hi.dot(xij);
            Scalar dWi = mKernel.gradValue(etai.magnitude(), detHi);
            Vector etaiHat = etai.unitVector();
            Vector gradWi = dWi * Hi.dot(etaiHat);

            Vector etaj = Hj.dot(xij);
            Scalar dWj = mKernel.gradValue(etaj.magnitude(), detHj);
            Vector etajHat = etaj.unitVector();
            Vector gradWj = dWj * Hj.dot(etajHat);
            Vector gradWij = 0.5 * (gradWi + gradWj);

            // Compute the Maxwell stress tensor.
            Scalar Bj2 = Bj.magnitude2();
            SymTensor Sj = (1.0/mMu0) * (Bj.selfdyad() - 
                                         0.5 * Bj2 * SymTensor(1., 0., 0.,
                                                               0., 1., 0.,
                                                               0., 0., 1.)) - mS0;

            // Assemble the contributions from the MHD equations.
            Vector aij;
//            Vector aij = mj * (Si.dot(gradWi)/(Omegai*rhoi*rhoi) +  // Non-symmetrized form
//                               Sj.dot(gradWj)/(Omegaj*rhoj*rhoj));
            if (mStabilizer == MorrisStabilizer)
              aij = -mj * (0.5*Bi.magnitude2()/(mMu0*rhoi*rhoi) + 
                           0.5*Bj.magnitude2()/(mMu0*rhoj*rhoj)) * gradWij +
                    mj/(mMu0*rhoi*rhoj) * (Bj.dyad(Bj) - Bi.dyad(Bi)).dot(gradWij);
            else
              aij = mj * (Si/(Omegai*rhoi*rhoi) + Sj/(Omegaj*rhoj*rhoj)) * gradWij; // Symmetrized form
            if (mStabilizer == B0rveStabilizer)
              aij -= Bi * mj * (Bj/(rhoj*rhoj) + Bi/(rhoi*rhoi)).dot(gradWi);
            dVdti += aij;
if (aij.magnitude() > maxAcc)
   maxAcc = aij.magnitude();
            dBdti -= mj/(Omegai*rhoi) * (Vij * Bi.dot(gradWi) - 
                                         Bi * Vij.dot(gradWi));

            // If we're using the compatible thermal energy discretization, 
            // feed it the pairwise accelerations.
            if (mSpecificThermalEnergyUpdate == computeCompatibleSpecificThermalEnergy)
            {
               pairAcci.push_back(aij);
               pairAccj.push_back(-mi*aij/mj);
            }

            // Total energy equation.
            if (mTotalSpecificEnergyUpdate == integrateTotalSpecificEnergy)
            {
               DeDt += mj * ((Si.dot(Vj) + Pi*Vj).dot(gradWi)/(Omegai*rhoi*rhoi) + 
                             (Sj.dot(Vi) + Pj*Vi).dot(gradWj)/(Omegaj*rhoj*rhoj));
            }
   
            // Compute the contribution to the gradient of the hyperbolic 
            // cleaning potential.
            gradPsi -= mj * (psii - psij) * gradWi / (Omegai * rhoi);
          }
        }
      } 

      if (mDivBCleaner == hyperbolicCleaner)
      {
         // Compute the contribution to the time derivative of the 
         // hyperbolic divergence cleaning potential.
         double ui = (*u[iNodeList])[iNode];
         double hi = pow(detHi, 1.0/3.0);
         double divBi = (*divB[iNodeList])[iNode];
         double gammai = condFluidPtr->equationOfState().gamma(rhoi, ui);
         assert(gammai*Pi/rhoi >= 0.0);
         double ch = sqrt(gammai*Pi/rhoi + 0.5*Bi2/(mMu0*rhoi));
         double sigma = 0.4; // Should be between 0.4 and 0.8.
         double tau = hi / (sigma*ch);
         if (DpsiDt.numFields() > 0)
         {
            CHECK(DpsiDt.numFields() == DvDt.numFields());
            (*DpsiDt[iNodeList])[iNode] = -ch*ch*divBi - psii/tau;
//printf("Dpsi/Dt[%zd][%d] = -%g**2*%g - %g/%g = %g\n", iNodeList, iNode, ch, divBi, psii, tau, -ch*ch*divBi - psii/tau);
         }

         // Add the hyperbolic divergence cleaning term in if necessary.
         dBdti -= gradPsi;
      }
    } 
  } // end for
//cout << "max acceleration: " << maxAcc << endl;

  // Now add resistive effects.
  if (isResistive) 
     mAddDiffusionDerivatives(dataBase, state, time, dt, derivs);

  // That's it.
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
MHD::
registerState(DataBase<Dim<3> >& dataBase,
              State<Dim<3> >& state) {

  typedef State<Dim<3> >::ScalarPolicyPointerType ScalarPolicyPointer;
  typedef State<Dim<3> >::VectorPolicyPointerType VectorPolicyPointer;

  // First, make sure that we've allocated storage for any fields we need.
  if (mDivBCleaner == hyperbolicCleaner)
  {
    if (mPsi.numFields() == 0)
    {
       for (DataBase<Dim<3> >::ConstFluidNodeListIterator 
            itr = dataBase.fluidNodeListBegin(); 
            itr != dataBase.fluidNodeListEnd(); ++itr)
       {
         mPsi.appendField(Field<Dim<3>, Dim<3>::Scalar>(MHDFieldNames::hyperbolicCleaning, **itr));
         mDpsiDt.appendField(Field<Dim<3>, Dim<3>::Scalar>(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::hyperbolicCleaning, **itr));
       }
    }
    CHECK(mPsi.numFields() == dataBase.numFluidNodeLists());
    CHECK(mDpsiDt.numFields() == dataBase.numFluidNodeLists());
  }

  // Traverse the set of fluid node lists in the database.
  size_t i = 0;
  for (DataBase<Dim<3> >::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++i) {

    // The Hydro should have set up the velocity and specific thermal energy
    // updates, so we don't have to.  We focus only on electromagnetic stuff, 
    // which requires that we use conducting fluids.
    ConductingFluidNodeList* condFluidPtr = 
      dynamic_cast<ConductingFluidNodeList*>(*itr);
    if (condFluidPtr == 0) continue;

    // Make some fake derivatives so that we can immediately invoke the 
    // update policies we need.
    vector<Physics<Dim<3> >*> dummyPackages;
    StateDerivatives<Dim<3> > derivs(dataBase, dummyPackages);

    // Magnetic field.
    VectorPolicyPointer BPolicy(new IncrementState<Dim<3>, Field<Dim<3>, Vector> >());
    state.registerField(condFluidPtr->magneticInduction(), BPolicy);

    // Magnetosonic speed.
    ScalarPolicyPointer CsPolicy(new MagnetosonicSpeedUpdatePolicy(mMu0));
    FieldList<Dim<3>, Scalar> soundSpeeds = state.scalarFields(HydroFieldNames::soundSpeed);
    FieldList<Dim<3>, Scalar>::iterator CsItr = soundSpeeds.fieldForNodeList(**itr);
    CHECK(CsItr < soundSpeeds.end());
    state.registerField(**CsItr, CsPolicy);
    const State<Dim<3> >::FieldKeyType CsKey(*itr, HydroFieldNames::soundSpeed);
    CsPolicy->update(CsKey, state, derivs, 1.0, 0.0, 1.0);

    // Resistivity.
    state.registerField(condFluidPtr->resistivity(), condFluidPtr->resistivityPolicy());

    // Current density.
    VectorPolicyPointer JPolicy(new CurrentDensityUpdatePolicy(mKernel, dataBase, mMu0));
    state.registerField(condFluidPtr->currentDensity(), JPolicy);

    // Divergence of the magnetic field, which is updated by the current 
    // density policy.
    ScalarPolicyPointer divBPolicy(new NonDynamicState<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >());
    state.registerField(condFluidPtr->magneticDivergence(), divBPolicy);

    // The hyperbolic divergence cleaning potential.
    if (mDivBCleaner == hyperbolicCleaner)
    {
      ScalarPolicyPointer psiPolicy(new IncrementState<Dim<3>, Field<Dim<3>, Scalar> >());
      state.registerField(*(mPsi[i]), psiPolicy);
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
MHD::
registerDerivatives(DataBase<Dim<3> >& dataBase,
                    StateDerivatives<Dim<3> >& derivs) {

  size_t i = 0;
  for (DataBase<Dim<3> >::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++i) {
    // Is this a conducting fluid?
    ConductingFluidNodeList* condFluidPtr = 
      dynamic_cast<ConductingFluidNodeList*>(*itr);
    if (condFluidPtr == 0) continue;

    derivs.registerField(condFluidPtr->DBDt());

    // The hyperbolic divergence cleaning potential.
    if (mDpsiDt.numFields() > 0)
       derivs.registerField(*(mDpsiDt[i]));
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
MHD::
applyGhostBoundaries(State<Dim<3> >& state,
                     StateDerivatives<Dim<3> >& derivs) {

  REQUIRE(this->valid());

  // Apply boundary conditions to the pertinent fields, assuming the 
  // Hydro does the same.
  FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);
  FieldList<Dim<3>, Scalar> divB = state.scalarFields(MHDFieldNames::magneticDivergence);
  FieldList<Dim<3>, Scalar> R = state.scalarFields(MHDFieldNames::resistivity);
  FieldList<Dim<3>, Vector> J = state.vectorFields(MHDFieldNames::currentDensity);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(B);
    (*boundaryItr)->applyFieldListGhostBoundary(divB);
    (*boundaryItr)->applyFieldListGhostBoundary(R);
    (*boundaryItr)->applyFieldListGhostBoundary(J);
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
MHD::
enforceBoundaries(State<Dim<3> >& state,
                  StateDerivatives<Dim<3> >& derivs) {

  REQUIRE(this->valid());

  // Enforce boundary conditions on the fields, assuming that the Hydro
  // has done its job.
  FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);
  FieldList<Dim<3>, Scalar> divB = state.scalarFields(MHDFieldNames::magneticDivergence);
  FieldList<Dim<3>, Scalar> R = state.scalarFields(MHDFieldNames::resistivity);
  FieldList<Dim<3>, Vector> J = state.vectorFields(MHDFieldNames::currentDensity);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(B);
    (*boundaryItr)->enforceFieldListBoundary(divB);
    (*boundaryItr)->enforceFieldListBoundary(R);
    (*boundaryItr)->enforceFieldListBoundary(J);
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void 
MHD::
initialize(const Scalar& time, 
           const Scalar& dt,
           const DataBase<Dim<3> >& db, 
           State<Dim<3> >& state,
           StateDerivatives<Dim<3> >& derivs)
{

  // If we are taking our first step, we need to compute the current density 
  // and the magnetic divergence, and apply our divergence cleaning algorithm
  // if necessary.
  if (mFirstStep)
  {
     // Make sure the current density and grad h are computed for each fluid.
     const FieldList<Dim<3>, Vector> J = state.vectorFields(MHDFieldNames::currentDensity);
     int numNodeLists = J.numFields();
     size_t i = 0;
     for (DataBase<Dim<3> >::ConstFluidNodeListIterator itr = db.fluidNodeListBegin();
          itr != db.fluidNodeListEnd(); ++itr, ++i) {
        ConductingFluidNodeList* condFluidPtr = 
           dynamic_cast<ConductingFluidNodeList*>(*itr);
        if (condFluidPtr == 0) continue;

        // J
        typedef State<Dim<3> >::VectorPolicyPointerType VectorPolicyPointer;
        VectorPolicyPointer JPolicy = state.policyForField(condFluidPtr->currentDensity());
        const State<Dim<3> >::FieldKeyType JKey(*itr, MHDFieldNames::currentDensity);
        JPolicy->update(JKey, state, derivs, 1.0, time, dt);

        // Omega
        typedef State<Dim<3> >::ScalarPolicyPointerType ScalarPolicyPointer;
        ScalarPolicyPointer OmegaPolicy = state.policyForField(condFluidPtr->omegaGradh());
        CHECK(OmegaPolicy.get() != 0);
        if (dynamic_cast<OmegaGradhPolicy<Dim<3> >*>(OmegaPolicy.get()) != 0)
        {
           const State<Dim<3> >::FieldKeyType OmegaKey(*itr, HydroFieldNames::omegaGradh);
           OmegaPolicy->update(OmegaKey, state, derivs, 1.0, time, dt);
        }
        else
        {
           // WTF?  I suppose I'll set Omega to 1 if I can't find the right 
           // update policy.
           condFluidPtr->omegaGradh() = 1.0;
        }
     }

     // Compute the divergence of B.
     mComputeDivB(db, state);

     // Clean the divergence.
     mCleanDivB(db, state);

     // Recompute div B.
     mComputeDivB(db, state);

     // FIXME: Recompute J?

     // We won't take our first step again.
     mFirstStep = false;
  }

  // If we are using the maximum stress stabilizer, compute the maximum 
  // stress to be subtracted from the stress terms in the induction equation.
  if (mStabilizer == maxStressStabilizer)
  {
    Scalar maxS = 0.0;

    // Get some field data.
    const FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);
    const FieldList<Dim<3>, Scalar> P = state.scalarFields(HydroFieldNames::pressure);

    // Now traverse the nodelists and recompute the maximum stress.
    int numNodeLists = B.numFields();
    for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
      // Is this a conducting fluid?  If not, we can skip it, since there's no 
      // magnetic field otherwise.
      const NodeList<Dim<3> >& nodeList = B[iNodeList]->nodeList(); 
      const ConductingFluidNodeList* condFluidPtr = 
        dynamic_cast<const ConductingFluidNodeList*>(&nodeList);
      if (condFluidPtr == 0) continue;
      for (int iNode = 0; iNode < nodeList.numInternalNodes(); ++iNode) {
        const Vector& Bi = (*B[iNodeList])[iNode];
        Scalar Pi = (*P[iNodeList])[iNode];
        Scalar Si = 0.5*Bi.magnitude2()/mMu0 - Pi;
        if (maxS < Si)
           maxS = Si;
      }
    }

    // Set the maximum stress.
    mS0 = SymTensor(maxS, maxS, maxS,
                    maxS, maxS, maxS,
                    maxS, maxS, maxS);
  }

  // Otherwise, if we are using the external field stabilizer, use the given 
  // value of the external field to determine the stress that will be 
  // subtracted.
  else if (mStabilizer == extFieldStabilizer)
  {
     mS0 = 1.0/mMu0 * mBext.selfdyad();
  }
  else
  {
     // Reset the stress attenuator to zero.
     mS0 = SymTensor();
  }

  // Compute the non-zero structure for the diffusion matrix unless all of 
  // our fluids are perfect conductors.
  if (hasResistivity(db, state))
  {
     mComputeMatrixStructure(db, state);
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void 
MHD::mComputeDivB(const DataBase<Dim<3> >& dataBase, 
                  State<Dim<3> >& state) const
{
   // Access to pertinent fields in the database.
   const FieldList<Dim<3>, Scalar> m = state.scalarFields(HydroFieldNames::mass);
   const FieldList<Dim<3>, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
   const FieldList<Dim<3>, Vector> x = state.vectorFields(HydroFieldNames::position);
   const FieldList<Dim<3>, Vector> v = state.vectorFields(HydroFieldNames::velocity);
   const FieldList<Dim<3>, Vector> J = state.vectorFields(MHDFieldNames::currentDensity);
   const FieldList<Dim<3>, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
   const FieldList<Dim<3>, Scalar> Omega = state.scalarFields(HydroFieldNames::omegaGradh);
   const FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);
   FieldList<Dim<3>, Scalar> divB = state.scalarFields(MHDFieldNames::magneticDivergence);

   // Reset the diagnostics for div B.
   Vector xMaxDivB;
   mMaxDivB = -FLT_MAX;
   mMinDivB = FLT_MAX;
   mAvgDivB = 0.0;
   size_t N = 0;

   // Compute!
   int numNodeLists = B.numFields();
   const ConnectivityMap<Dim<3> >& connectivityMap = dataBase.connectivityMap();
   for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
      const NodeList<Dim<3> >& nodeList = J[iNodeList]->nodeList(); 
      for (int iNode = 0; iNode < nodeList.numInternalNodes(); ++iNode) {
         // Get data for the ith node.
         Scalar mi = (*m[iNodeList])[iNode];
         Scalar rhoi = (*rho[iNodeList])[iNode];
         const Vector& Bi = (*B[iNodeList])[iNode];
         const SymTensor& Hi = (*H[iNodeList])[iNode];
         Scalar detHi = Hi.Determinant();
         Scalar Omegai = (*Omega[iNodeList])[iNode];
         Scalar divBi = 0.;

         // Iterate over the neighboring NodeLists.
         const vector< vector<int> >& fullConnectivity = 
            connectivityMap.connectivityForNode(&nodeList, iNode);
         for (int jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
            const vector<int>& connectivity = fullConnectivity[jNodeList];
            if (connectivity.size() > 0) {
               // Loop over the neighbors.
               for (vector<int>::const_iterator jItr = connectivity.begin();
                     jItr != connectivity.end();
                     ++jItr) {
                  int jNode = *jItr;
                  Vector xij = (*x[iNodeList])[iNode] - (*x[jNodeList])[jNode];
                  CHECK(xij.magnitude2() != 0.0);
                  Scalar mj = (*m[jNodeList])[jNode];
                  const Vector& Bj = (*B[jNodeList])[jNode];
                  Vector etai = Hi.dot(xij);
                  Scalar dWi = mKernel.gradValue(etai.magnitude(), detHi);
                  Vector etaiHat = etai.unitVector();
                  Vector gradWi = dWi * Hi.dot(etaiHat);
                  divBi -= mj/(rhoi*Omegai) * (Bi - Bj).dot(gradWi);
               }
            }
         }
         (*divB[iNodeList])[iNode] = divBi;

         // Update the local diagnostics for div B.
         if (divBi > mMaxDivB)
         {
            mMaxDivB = divBi;
            xMaxDivB = (*x[iNodeList])[iNode];
         }
         if (divBi < mMinDivB)
            mMinDivB = divBi;
         mAvgDivB += divBi;
      }
      N += nodeList.numInternalNodes();
   } // end for

   // Now share the diagnostics across all processors.
   double globalVal, localVals[2], globalVals[2];
   MPI_Allreduce(&mMaxDivB, &globalVal, 1, MPI_DOUBLE, MPI_MAX, Communicator::communicator());
   mMaxDivB = globalVal;
   MPI_Allreduce(&mMinDivB, &globalVal, 1, MPI_DOUBLE, MPI_MIN, Communicator::communicator());
   mMinDivB = globalVal;
   localVals[0] = mAvgDivB, localVals[1] = 1.0 * N;
   MPI_Allreduce(localVals, globalVals, 2, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
   mAvgDivB = globalVals[0]/globalVals[1];

#if 0
int rank;
MPI_Comm_rank(Communicator::communicator(), &rank);
if (rank == 0)
{
   printf("Max div B = %g, min div B = %g, avg div B = %g\n", mMaxDivB, mMinDivB, mAvgDivB);
   printf("(max div B at x = (%g, %g, %g)\n", xMaxDivB.x(), xMaxDivB.y(), xMaxDivB.z());
}
#endif
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void 
MHD::mCleanDivB(const DataBase<Dim<3> >& dataBase, 
                State<Dim<3> >& state) const
{
   // Access to pertinent fields in the database.
   const FieldList<Dim<3>, Scalar> m = state.scalarFields(HydroFieldNames::mass);
   const FieldList<Dim<3>, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
   const FieldList<Dim<3>, Vector> x = state.vectorFields(HydroFieldNames::position);
   const FieldList<Dim<3>, Vector> v = state.vectorFields(HydroFieldNames::velocity);
   const FieldList<Dim<3>, Vector> J = state.vectorFields(MHDFieldNames::currentDensity);
   const FieldList<Dim<3>, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
   const FieldList<Dim<3>, Scalar> Omega = state.scalarFields(HydroFieldNames::omegaGradh);
   const FieldList<Dim<3>, Scalar> divB = state.scalarFields(MHDFieldNames::magneticDivergence);

   // Fields we'll be modifying.
   FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);

   const ConnectivityMap<Dim<3> >& connectivityMap = dataBase.connectivityMap();
   int numNodeLists = B.numFields();

   if ((mDivBCleaner == GreensFnProjCleaner) ||
       (mDivBCleaner == BiotSavartProjCleaner))
   {
      // Compute the corrected B.
      for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
         const NodeList<Dim<3> >& nodeList = B[iNodeList]->nodeList(); 

         // Is this a conducting fluid?  If not, we can skip it, since there's no 
         // magnetic field otherwise.
         const ConductingFluidNodeList* condFluidPtr = 
            dynamic_cast<const ConductingFluidNodeList*>(&nodeList);
         if (condFluidPtr == 0) continue;

         // Now compute those derivatives.
         for (int iNode = 0; iNode < nodeList.numInternalNodes(); ++iNode) {
            Vector& Bi = (*B[iNodeList])[iNode];
            const Vector& xi = (*x[iNodeList])[iNode];

            // Here are various versions of to-be-corrected magnetic inductions.
            Vector Bg = Bi;      // Greens function projection
            Vector Bbs;          // Biot-Savart projection

            // Loop over the neighboring nodes.
            const vector< vector<int> >& fullConnectivity = 
               connectivityMap.connectivityForNode(&nodeList, iNode);
            // Iterate over the NodeLists.
            for (int jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
               const vector<int>& connectivity = fullConnectivity[jNodeList];
               if (connectivity.size() > 0) {
                  // Loop over the neighbors.
                  for (vector<int>::const_iterator jItr = connectivity.begin();
                        jItr != connectivity.end();
                        ++jItr) {
                     int jNode = *jItr;
                     CHECK(iNode != jNode);

                     // Nodal data.
                     Scalar mj = (*m[jNodeList])[jNode];
                     Scalar rhoj = (*rho[jNodeList])[jNode];
                     Vector xij = xi - (*x[jNodeList])[jNode];
                     double xijMag = xij.magnitude();
                     CHECK(xijMag != 0.0);
                     const Vector& Jj = (*J[jNodeList])[jNode];
                     Vector curlBj = mMu0 * Jj;
                     double divBj = (*divB[jNodeList])[jNode];

                     // Green's function projection.
                     Bg += mj * divBj * xij / (4*M_PI*rhoj*xijMag*xijMag*xijMag);
//printf("Bg[%d] += %g * %g * (%g, %g, %g)/(4 * pi * %g * %g**3) -> (%g, %g, %g)\n", iNode, mj, divBj, xij.x(), xij.y(), xij.z(), rhoj, xijMag, Bg.x(), Bg.y(), Bg.z());

                     // Biot-Savart projection.
                     Bbs -= mj * curlBj.cross(xij)/(4*M_PI*rhoj*xijMag*xijMag*xijMag);
//printf("Bbs[%d] -= %g * (%g, %g, %g) x (%g, %g, %g)/(4 * pi * %g * %g**3) -> (%g, %g, %g)\n", iNode, mj, curlBj.x(), curlBj.y(), curlBj.z(), xij.x(), xij.y(), xij.z(), rhoj, xijMag, Bbs.x(), Bbs.y(), Bbs.z());
                  }
               }
            } 

            // Replace the magnetic induction with its "cleaned" value.
            if (mDivBCleaner == GreensFnProjCleaner)
               Bi = Bg;
            else if (mDivBCleaner == BiotSavartProjCleaner)
               Bi = Bbs;
         } 
      } // end for
   }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void 
MHD::postStateUpdate(const DataBase<Dim<3> >& dataBase, 
                     State<Dim<3> >& state,
                     const StateDerivatives<Dim<3> >& derivatives) const
{
   // Compute the divergence of B. 
   mComputeDivB(dataBase, state);

   // Clean the divergence and recompute it.
   mCleanDivB(dataBase, state);
   mComputeDivB(dataBase, state);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void 
MHD::finalize(const Scalar time, 
              const Scalar dt,
              DataBase<Dim<3> >& dataBase, 
              State<Dim<3> >& state,
              StateDerivatives<Dim<3> >& derivs)
{
  // Tally the magnetic field energy.
  mMagneticEnergy = 0.0;
  const FieldList<Dim<3>, Vector> B = state.vectorFields(MHDFieldNames::magneticInduction);
  const FieldList<Dim<3>, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dim<3>, Scalar> m = state.scalarFields(HydroFieldNames::mass);

  // Now traverse the nodelists and recompute the maximum stress.
  int numNodeLists = B.numFields();
  for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
    // Is this a conducting fluid?  If not, we can skip it, since there's no 
    // magnetic field otherwise.
    const NodeList<Dim<3> >& nodeList = B[iNodeList]->nodeList(); 
    const ConductingFluidNodeList* condFluidPtr = 
      dynamic_cast<const ConductingFluidNodeList*>(&nodeList);
    if (condFluidPtr == 0) continue;
    for (int iNode = 0; iNode < nodeList.numInternalNodes(); ++iNode) {
      Scalar mi = (*m[iNodeList])[iNode];
      Scalar rhoi = (*rho[iNodeList])[iNode];
      const Vector& Bi = (*B[iNodeList])[iNode];
      mMagneticEnergy += 0.5 * mi * Bi.magnitude2()/(rhoi*mMu0);
    }
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
MHD::TimeStepType
MHD::
dt(const DataBase<Dim<3> >& dataBase, 
   const State<Dim<3> >& state,
   const StateDerivatives<Dim<3> >& derivs,
   const Scalar currentTime) const
{
  return TimeStepType(FLT_MAX, "No additional timestep constraints.");
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
MHD::Scalar 
MHD::
extraEnergy() const
{
  return mMagneticEnergy;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
boost::python::handle<PyObject>
MHD::
diffusionMatrix() const
{
  return boost::python::handle<PyObject>(mDiffMatrix);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
boost::python::handle<PyObject>
MHD::
diffusionRHS() const
{
  return boost::python::handle<PyObject>(mDiffRHS);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
MHD::
implicitness(double alpha)
{
  REQUIRE((alpha >= 0.0) && (alpha <= 1.0));
  mImplicitness = alpha;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
bool 
MHD::
valid() const
{
#if 0
  // Make sure that we have only parallel and periodic boundary 
  // conditions for now.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    if ((dynamic_cast<DistributedBoundary<Dim<3> >*>(*boundaryItr) == 0) and 
        (dynamic_cast<PeriodicBoundary<Dim<3> >*>(*boundaryItr) == 0))
    {
       return false;
    }
  }
#endif
  return true;
}

//------------------------------------------------------------------------------
} // end namespace Spheral
