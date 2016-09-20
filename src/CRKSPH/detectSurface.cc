//---------------------------------Spheral++------------------------------------
// Detect surface particles leveraging the zeroth and first moments
//------------------------------------------------------------------------------
#include <stdio.h>

#include "detectSurface.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
    namespace CRKSPH {
        using namespace std;
        using std::min;
        using std::max;
        using std::abs;
        
        using FieldSpace::Field;
        using FieldSpace::FieldList;
        using NeighborSpace::ConnectivityMap;
        using KernelSpace::TableKernel;
        using NodeSpace::NodeList;
        
        template<typename Dimension>
        void
        detectSurface(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                      const double detectThreshold,
                      const double detectRange,
                      const double sweepAngle,
                      FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& surfNorm) {
            // Pre-conditions.
            const size_t numNodeLists = m0.size();
            REQUIRE(m1.size() == numNodeLists);
            REQUIRE(surfNorm.size() == numNodeLists);
            
            typedef typename Dimension::Scalar Scalar;
            typedef typename Dimension::Vector Vector;
            typedef typename Dimension::Tensor Tensor;
            
            const Scalar nPerh = m0.nodeListPtrs()[0]->nodesPerSmoothingScale();
            
            // Zero out surfNorm to get fresh maps at each step
            surfNorm = 0;
            
            // Walk the FluidNodeLists.
            for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
                const NodeList<Dimension>& nodeList = m0[nodeListi]->nodeList();
                const int firstGhostNodei = nodeList.firstGhostNode();
                
                // Iterate over nodes in nodeListi
                for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
                     iItr != connectivityMap.end(nodeListi);
                     ++iItr) {
                    const int i = *iItr;
                    
                    // Get the state for node i.
                    const Vector& ri    = position(nodeListi, i);
                    const SymTensor& Hi = H(nodeListi, i);
                    const Scalar m0i    = m0(nodeListi, i);
                    const Vector& m1i   = m1(nodeListi, i);
                    const Vector m1ih   = m1i.unitVector();
                    bool particleDetected = 0;
                    
                    if (m0i < detectThreshold) {
                        // Get neighbors
                        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
                        CHECK(fullConnectivity.size() == numNodeLists);
                        // Loop over them
                        for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
                            if (particleDetected) break;
                            const vector<int>& connectivity = fullConnectivity[nodeListj];
                            const int firstGhostNodej = m0[nodeListj]->nodeList().firstGhostNode();
                            for (vector<int>::const_iterator jItr = connectivity.begin();
                                 jItr != connectivity.end();
                                 ++jItr) {
                                const int j = *jItr;
                                
                                // Get the state for node j.
                                const Vector& rj    = position(nodeListj, j);
                                const SymTensor& Hj = H(nodeListj, j);
                                const Vector rji    = rj - ri;  // pointing away from i as does m1
                                const Vector rjih   = rji.unitVector();
                                
                                // Check range
                                
                                // Check angle
                                const cangle = m1ih.dot(rjih);  // because they're both unit vectors!
                                if ((cangle >= cos(sweepAngle)) && (rji.magnitude() < detectRange)) {
                                    particleDetected = 1;
                                    break;
                                }// needs to be scaled to eta space
                            }
                        }
                    }
                    
                    if !(particleDetected) surfNorm(nodeListi, i) = 1;
                }
            }
        }
    }
}