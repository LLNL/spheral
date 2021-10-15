from Spheral3d import *
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
for p in [Vector(1, 0, 0), Vector(1, 1, 0), Vector(0, 1, 0), Vector(0, 0, 0),
          Vector(1, 0, 1), Vector(1, 1, 1), Vector(0, 1, 1), Vector(0, 0, 1),
          Vector(0.5, 1.5, 1)]:
    nodePositions.append(p)

# Edges.
edgeNodes = list2vec([(0, 1), (1, 2), (2, 3), (3, 0),
                      (0, 4), (1, 5), (2, 6), (3, 7),
                      (4, 5), (5, 6), (6, 7), (7, 4),
                      (1, 8), (2, 8), (5, 8), (6, 8)])

# Faces.
faceEdges = list2vec([(0, 1, 2, 3), (0, 5, 8, 4), (2, 6,10, 7), (3, 7,11, 4),
                      (8, 9, 10, 11), (5, 12, 14), (1, 13, 12), (13, 15, 6), (9, 14, 15)])
                      

# Zones.
zoneFaces = list2vec([[0, 1, 2, 3, 4, 5, 6, 7, 8]])

mesh = PolyhedralMesh(nodePositions, edgeNodes, faceEdges, zoneFaces)
siloMeshDump("polyhedral_example", mesh, label = "Test polyhedral mesh.")
