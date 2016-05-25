from Spheral2d import *
from siloMeshDump import *

def list2vec(stuff):
    result = vector_of_vector_of_unsigned()
    for x in stuff:
        result.append(vector_of_unsigned())
        for y in x:
            result[-1].append(y)
    return result

# Nodes.
nodePositions = vector_of_Vector()
for p in [Vector(0, 0),
          Vector(1, 0),
          Vector(0, 1),
          Vector(1, 1)]:
    nodePositions.append(p)

# Edges.
edgeNodes = list2vec([[0, 1],
                      [1, 2],
                      [2, 3],
                      [3, 0]])

# Faces.
faceEdges = list2vec([[0], [1], [2], [3]])

# Zones.
zoneFaces = list2vec([[0, 1, 2, 3]])

mesh = PolygonalMesh(nodePositions, edgeNodes, faceEdges, zoneFaces)
siloMeshDump("polygonal_example", mesh, label = "Test polygonal mesh.")
