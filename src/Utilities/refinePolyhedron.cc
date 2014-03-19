//---------------------------------Spheral++----------------------------------//
// refinePolyhedron
//
// Refine/smooth a polyhedron.
//
// We use the Pixar OpenSubdiv package for our work here.  This implementation
// is based on an example I found in their source:  examples/tessellateObjFile.
//
// Created by JMO, Tue Mar 18 15:23:04 PDT 2014
//----------------------------------------------------------------------------//
#include "refinePolyhedron.hh"

// We base this on Pixar's opensubdiv package.
#include "opensubdiv/osdutil/adaptiveEvaluator.h"
#include "opensubdiv/osdutil/uniformEvaluator.h"
#include "opensubdiv/osdutil/topology.h"
#include "opensubdiv/osd/error.h"

#include <vector>
#include <algorithm>

#include "Geometry/Dimension.hh"

namespace Spheral {

using namespace std;

GeomPolyhedron refinePolyhedron(const GeomPolyhedron& poly0,
                                const unsigned numLevels) {

  typedef Dim<3>::Vector Vector;
  typedef GeomPolyhedron::Facet Facet;

  // Start by reading out the Spheral polyhedral data to arrays we're going
  // to shove into opensubdiv.
  const int numVertices0 = poly0.vertices().size();
  const int numFaces0 = poly0.facets().size();
  vector<float> positions0(3*numVertices0);
  vector<int> faceVertices0, numVertsPerFace0(numFaces0);
  {
    const vector<Vector>& verts0 = poly0.vertices();
    const vector<Facet>& facets0 = poly0.facets();
    for (unsigned i = 0; numVertices0; ++i) {
      positions0[3*i  ] = verts0[i].x();
      positions0[3*i+1] = verts0[i].y();
      positions0[3*i+2] = verts0[i].z();
    }
    for (unsigned i = 0; i != numFaces0; ++i) {
      const vector<unsigned>& ipoints = facets0[i].ipoints();
      numVertsPerFace0[i] = ipoints.size();
      for (unsigned j = 0; j != ipoints.size(); ++j) {
        faceVertices0.push_back(ipoints[j]);
      }
    }
  }
    
  // Initialize the OSD topology.
  PxOsdUtilSubdivTopology topology;
  bool ok;
  std::string errorMessage;
  ok = topology.Initialize(numVertices0, 
                           &numVertsPerFace0[0], numVertsPerFace0.size(),
                           &faceVertices0[0], faceVertices0.size(),
                           numLevels,
                           &errorMessage);
  VERIFY2(ok, errorMessage);

  // Create a uniform refinement thingus.
  PxOsdUtilUniformEvaluator uniformEvaluator;
  ok = uniformEvaluator.Initialize(topology, &errorMessage);
  VERIFY2(ok, errorMessage);

  // Give the evaluator the intial vertex positions to start with.
  uniformEvaluator.SetCoarsePositions(positions0, &errorMessage);

  // Do some refining (the "1" here means single thread).
  ok = uniformEvaluator.Refine(1, &errorMessage);
  VERIFY2(ok, errorMessage);

  // Extract the refined topology and positions.
  PxOsdUtilSubdivTopology topology1;
  const float *positions1 = NULL;
  ok = uniformEvaluator.GetRefinedTopology(&topology1, &positions1, &errorMessage);
  VERIFY2(ok, errorMessage);

  const int numVertices1 = topology1.numVertices;
}

}
