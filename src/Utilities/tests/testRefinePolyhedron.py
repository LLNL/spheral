from Spheral3d import *

# Create a simple unit box.
verts0 = vector_of_Vector()
verts0.append(Vector(0,0,0))
verts0.append(Vector(1,0,0))
verts0.append(Vector(1,1,0))
verts0.append(Vector(0,1,0))
verts0.append(Vector(0,0,1))
verts0.append(Vector(1,0,1))
verts0.append(Vector(1,1,1))
verts0.append(Vector(0,1,1))
poly0 = Polyhedron(verts0)

# Now refine it.
poly1 = refinePolyhedron(poly0, 1)


