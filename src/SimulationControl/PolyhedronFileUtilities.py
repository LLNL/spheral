#-------------------------------------------------------------------------------
# A set of methods to help reading and writing Spheral Polyhedra to/from files.
#-------------------------------------------------------------------------------
from Spheral3d import Vector, Tensor, SymTensor, Polyhedron, \
    vector_of_Vector, vector_of_unsigned, vector_of_vector_of_unsigned

#-------------------------------------------------------------------------------
# Read an OBJ (vertex-facet labeled) shape file to create a polyhedron.
# This is the format used by the NASA radar asteroid shape models:
#  http://echo.jpl.nasa.gov/asteroids/shapes/shapes.html
# NOTE!!  These things seem to count from 1 rather than 0 for indices!
#-------------------------------------------------------------------------------
def readPolyhedronOBJ(filename):
    f = open(filename, "r")
    verts = vector_of_Vector()
    facets = vector_of_vector_of_unsigned()
    for line in f:
        stuff = line.split()
        if stuff[0] == "v":
            assert len(stuff) == 4
            verts.append(Vector(float(stuff[1]), float(stuff[2]), float(stuff[3])))
        elif stuff[0] == "f":
            assert len(stuff) >= 4
            facets.append(vector_of_unsigned())
            for x in stuff[1:]:
                facets[-1].append(int(x) - 1)
    f.close()
    nverts = len(verts)
    for i in xrange(len(facets)):
        for j in xrange(len(facets[i])):
            assert facets[i][j] < nverts
    poly = Polyhedron(verts, facets)
    return poly

#-------------------------------------------------------------------------------
# Write an OBJ (vertex-facet labeled) shape file from a polyhedron.
# NOTE!!  These things seem to count from 1 rather than 0 for indices!
#-------------------------------------------------------------------------------
def writePolyhedronOBJ(poly, filename, forceTriangles=False):
    f = open(filename, "w")
    verts = poly.vertices()
    facets = poly.facets()
    for v in verts:
        f.write("v %g %g %g\n" % (v.x, v.y, v.z))
    for facet in facets:
        ipoints = facet.ipoints
        f.write("f")
        if forceTriangles and len(ipoints > 3):
            for j in xrange(1, len(ipoints)-1):
                f.write(" %i %i %i\n" % (0, j, j+1))
        else:
            for i in ipoints:
                f.write(" %i" % (i + 1))
        f.write("\n")
    f.close()
    return

#-------------------------------------------------------------------------------
# Read an OFF shape file to create a polyhedron.  Used by geomview.
# http://people.sc.fsu.edu/~jburkardt/data/off/off.html
#-------------------------------------------------------------------------------
def readPolyhedronOFF(filename):
    f = open(filename, "r")
    line = f.readline()
    assert line[:3].lower() == "off"
    while line[0] == "#":
        line = f.readline()
    stuff = line.split()
    assert len(stuff) == 3
    nv, nf, ne = int(stuff[0]), int(stuff[1]), int(stuff[2])

    # Read the vertex positions.
    verts = vector_of_Vector()
    for i in xrange(nv):
        line = f.readline()
        stuff = line.split()
        assert len(stuff) == 3
        verts.append(Vector(float(stuff[0]), float(stuff[1]), float(stuff[2])))
    
    # Read the face vertex indices.
    facets = vector_of_vector_of_unsigned()
    for i in xrange(nf):
        line = f.readline()
        stuff = line.split()
        assert len(stuff) >= 4
        n = int(stuff[0])
        facets.push_back(vector_of_unsigned(n))
        for j in xrange(n):
            k = int(stuff[j+1])
            assert k < nv
            facets[-1][j] = k
    
    # Build the polyhedron and we're done.
    poly = Polyhedron(verts, facets)
    return poly


#-------------------------------------------------------------------------------
# Write an OFF shape file from a polyhedron.
# We fudge this a bit and don't write the edges number out.  Hopefully that 
# doesn't bust anything!
#-------------------------------------------------------------------------------
def writePolyhedronOFF(poly, filename):
    f = open(filename, "w")
    verts = poly.vertices()
    facets = poly.facets()
    f.write("OFF\n")
    f.write("%i %i %i\n" % (verts.size(), facets.size(), 0))

    # Write the vertex coordinates.
    for v in verts:
        f.write("%g %g %g\n" % (v.x, v.y, v.z))

    # Write the facet vertex indices.
    for facet in facets:
        ipoints = facet.ipoints
        f.write("%i" % ipoints.size())
        for i in ipoints:
            f.write(" %i" % i)
        f.write("\n")

    # That's it.
    f.close()
    return
