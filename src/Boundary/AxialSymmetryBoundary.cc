#include "AxialSymmetryBoundary.hh"

#include "Field/NodeIterators.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {


//------------------------------------------------------------------------------
AxialSymmetryBoundary::
AxialSymmetryBoundary(TableKernel<Dim<3> >* kernel):
   Boundary<Dim<3> >(),
   mKernel(kernel),
   mR() {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
AxialSymmetryBoundary::
~AxialSymmetryBoundary() {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
setGhostNodes(NodeList<Dim<3> >& nodeList) {
   // Add this NodeList, creating space for control & ghost nodes.
   addNodeList(nodeList);

   // Define control nodes (and here, every internal node is a control node).
   Boundary<Dim<3> >::BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
   vector<int>& controlNodes = boundaryNodes.controlNodes;
   vector<int>& ghostNodes = boundaryNodes.ghostNodes;
   controlNodes.resize(0);
   ghostNodes.resize(0);

   // If we have no internal nodes, bug out.
   if (nodeList.numInternalNodes() == 0) return;

   // We want to make a cylinder of ghost nodes around the existing (colinear) 
   // internal nodes.  We characterize this cylinder by its total radius R, 
   // its radial spacing, dr, and its angular spacing dPhi.

   // First of all, R is the total radial extent of the kernel.  This is 
   // equal to sqrt(hx^2 + hy^2), where hx and hy are the respective x- and y-
   // extents of the H tensor.  Since we are working with a line of nodes 
   // along the (z) axis, however, it's likely that these extents are 
   // completely made-up.  So we will use the approximation that 
   // hx = hy = hz, since hz should have been set properly, and this will 
   // give use a uniform distribution of ghost nodes.

   // Let M = H*H.  Then hz = sqrt(cofactor(M, 3, 3))/ det H.  Then 
   // R = sqrt(2 * hz * hz).
   Dim<3>::SymTensor H = nodeList.Hfield()(0);
   double detH = H.Determinant();
   double M11 = H(0,0)*H(0,0) + H(0,1)*H(0,1) + H(0,2)*H(0,2);
   double M21 = H(0,0)*H(0,1) + H(0,1)*H(1,1) + H(0,2)*H(1,2);
   double M31 = H(0,0)*H(0,2) + H(0,1)*H(1,2) + H(0,2)*H(2,2);
   double M22 = H(0,1)*H(0,1) + H(1,1)*H(1,1) + H(1,2)*H(1,2);
   double M32 = H(0,1)*H(0,2) + H(1,1)*H(1,2) + H(1,2)*H(2,2);
   double M33 = H(0,2)*H(0,2) + H(1,2)*H(1,2) + H(2,2)*H(2,2);
   double Mcof33 = M11*M22 - M21*M21;
   double h_z = sqrt(Mcof33) / detH;
   double R = sqrt(2.0*h_z*h_z);
   mR[&nodeList] = R;

   // Now R = Nr*dr, where N is the number of neighbors in the radial extent of 
   // the kernel.  Thus, Nr = (# of neighbors per H) * (# of H per kernel extent).
   // Below, we add 1 to Nr to make sure that we're creating enough neighbors.
   // Note that dr = R / Nr.
   int Nr = static_cast<int>(nodeList.nodesPerSmoothingScale() * mKernel->kernelExtent()) + 1;
   double dr = R / Nr;

   // Now set up the ghost nodes for this boundary condition.
   // We create a bevy of nodes around each internal node.
   int numberOfGhostNodesToAdd = 0;
   for (InternalNodeIterator<Dim<3> > iter = nodeList.internalNodeBegin();
         iter != nodeList.internalNodeEnd(); ++iter)
   {
      // Add this internal node as a control node.
      controlNodes.push_back(iter.nodeID());

      // Find the position of the internal node.
      Dim<3>::Vector xInternal = nodeList.positions()(iter);

      // We'll progress along the radial direction, sweeping out the azimuth.
      for (int i = 1; i < Nr; ++i) 
      {
         // Compute our radial position.
         double r = i*dr;

         // We choose angular resolution by equating r * dPhi and dr.  This 
         // means that dPhi = (dr / r).
         double dPhi = dr / r;

         // Progress around the azimuth.  There are Nphi neighbors around a 
         // given azimuth.
         int Nphi = static_cast<int>(2*M_PI / dPhi);
         numberOfGhostNodesToAdd += Nphi;
      } // end for
   } // end for

   // Allocate space for the new ghost nodes we're going to create.
   int currentNumGhostNodes = nodeList.numGhostNodes();
   int firstNewGhostNode = nodeList.numNodes();
   nodeList.numGhostNodes(currentNumGhostNodes + numberOfGhostNodesToAdd);
   ghostNodes.resize(numberOfGhostNodesToAdd);

   // Denote the ghost nodes.
   for (int i = 0; i < ghostNodes.size(); ++i)
   {
      ghostNodes[i] = firstNewGhostNode + i;
   } // end for

   // Update the ghost nodes.
   updateGhostNodes(nodeList);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
updateGhostNodes(NodeList<Dim<3> >& nodeList) {
   // We are going to update all the positions of our ghost nodes.  We 
   // do this by iterating over the internal nodes and updating the 
   // ghost nodes that correspond to these internal nodes.

   // If there are no internal nodes, bug out.
   if (nodeList.numInternalNodes() == 0) return;

   BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
   vector<int>& controlNodes = boundaryNodes.controlNodes;
   vector<int>& ghostNodes = boundaryNodes.ghostNodes;

   // Now iterate over the control nodes and set the positions of the 
   // ghost nodes for each control node.
   vector<int>::const_iterator controlItr = controlNodes.begin();
   vector<int>::const_iterator ghostItr = ghostNodes.begin();
   double R = mR[&nodeList];
   int Nr = static_cast<int>(nodeList.nodesPerSmoothingScale() * mKernel->kernelExtent()) + 1;
   double dr = R / Nr;
   for (; controlItr < controlNodes.end(); ++controlItr)
   {
      // Find the position of the internal node.
      Dim<3>::Vector xInternal = nodeList.positions()(*controlItr);

      // We'll progress along the radial direction, sweeping out the azimuth.
      for (int i = 1; i < Nr; ++i) 
      {
         // Compute our radial position.
         double r = i*dr;

         // We choose angular resolution by equating r * dPhi and dr.  This 
         // means that dPhi = (dr / r).
         double dPhi = dr / r;

         // Progress around the azimuth.  There are Nphi neighbors around a 
         // given azimuth.
         int Nphi = static_cast<int>(2*M_PI / dPhi);
         for (int j = 0; j < Nphi; ++j)
         {
            // The z-component of the ghost node is the same as that for the 
            // internal node.  The x and y ones are given by the (r, phi) 
            // coordinates deduced from the above (i, j) pair.
            double phi = j*(2.0 * M_PI / Nphi);
            Dim<3>::Vector xGhost(xInternal(0) + r * cos(phi),
                                  xInternal(1) + r * sin(phi),
                                  xInternal(2));
            nodeList.positions()(*ghostItr) = xGhost;
            ++ghostItr;
         } // end for
      } // end for
   } // end for

  // Set the Hfield.
  Field<Dim<3>, SymTensor>& Hfield = nodeList.Hfield();
  this->applyGhostBoundary(Hfield);

//   // Let the neighbor know we've updated.
//   nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  Technically we flag all internal nodes as "violation" nodes,
// since we want to make sure they all meet our fake "1-D" approximation.
//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
setViolationNodes(NodeList<Dim<3> >& nodeList) {
  // Loop through all the internal nodes and add them as violation nodes.
  Boundary<Dim<3> >::BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  vector<int>& violationNodes = boundaryNodes.violationNodes;
  for (InternalNodeIterator<Dim<3> > iter = nodeList.internalNodeBegin();
       iter != nodeList.internalNodeEnd(); ++iter)
  {
     violationNodes.push_back(iter.nodeID());
  } // end for
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
updateViolationNodes(NodeList<Dim<3> >& nodeList) {
  // Zero out the x and y components of every violation node.  We can do this
  // by looping through the violation nodes and just setting the positions 
  // of the nodes accordingly.
  Boundary<Dim<3> >::BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  vector<int>& violationNodes = boundaryNodes.violationNodes;
  Field<Dim<3>, Vector>& x = nodeList.positions();
  for (size_t i = 0; i < violationNodes.size(); ++i)
  {
    x[violationNodes[i]].x(0.0);
    x[violationNodes[i]].y(0.0);
  } // end for
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template <typename Value>
void
AxialSymmetryBoundary::
mApplyGhostBoundary(Field<Dim<3>, Value>& field) const {
   // Copy the field value from each control node to its ghost nodes.
   NodeList<Dim<3> >& nodeList = const_cast<NodeList<Dim<3> >&>(field.nodeList());
   const vector<int>& controlNodes = this->controlNodes(nodeList);
   const vector<int>& ghostNodes = this->ghostNodes(nodeList);
   double R = (mR.find(&nodeList))->second;
   int Nr = static_cast<int>(nodeList.nodesPerSmoothingScale() * mKernel->kernelExtent()) + 1;
   double dr = R / Nr;
   int lastGhostUpdated = 0;
   for (size_t i = 0; i < controlNodes.size(); ++i)
   {
      // Figure out how many ghost nodes there are for this control node.
      // We'll progress along the radial direction, sweeping out the azimuth.
      int numGhosts = 0;
      for (int j = 0; j < Nr; ++j) 
      {
         double r = j*dr;
         double dPhi = dr / r;
         int Nphi = static_cast<int>(2.0*M_PI / dPhi);
         numGhosts += Nphi;
      } // end for

      // Now update each of the values for the ghost nodes based upon the 
      // values at the control node.
      CHECK(ghostNodes[lastGhostUpdated + numGhosts - 1] <= field.size());
      for (int j = 0; j < numGhosts; ++j)
      {
         field[ghostNodes[lastGhostUpdated + j]] = field[controlNodes[i]];  
      } // end for
      lastGhostUpdated += numGhosts;

   } // end for
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, int>& field) const {
   mApplyGhostBoundary(field);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Scalar>& field) const {
   mApplyGhostBoundary(field);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {
   mApplyGhostBoundary(field);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {
   mApplyGhostBoundary(field);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {
   mApplyGhostBoundary(field);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field) const {
   mApplyGhostBoundary(field);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, std::vector<Dim<3>::Scalar> >& field) const {
   mApplyGhostBoundary(field);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
enforceBoundary(Field<Dim<3>, int>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Scalar>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
AxialSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field) const {
}
//------------------------------------------------------------------------------

}
