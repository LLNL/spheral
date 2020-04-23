#-------------------------------------------------------------------------------
# A set of methods to help reading and writing Spheral Polyhedra to/from files.
#-------------------------------------------------------------------------------
from SpheralCompiledPackages import *
import numpy as np
import stl

#-------------------------------------------------------------------------------
# Read an OBJ (vertex-facet labeled) shape file to create a polyhedron.
# This is the format used by the NASA radar asteroid shape models:
#  http://echo.jpl.nasa.gov/asteroids/shapes/shapes.html
# NOTE!!  These things seem to count from 1 rather than 0 for indices!
#-------------------------------------------------------------------------------
def readPolyhedronOBJ(filename):
    f = open(filename, "r")
    verts = []
    facets = []
    for line in f:
        stuff = line.split()
        if stuff:
            if stuff[0] == "v":
                assert len(stuff) == 4
                verts.append(Vector3d(float(stuff[1]), float(stuff[2]), float(stuff[3])))
            elif stuff[0] == "f":
                assert len(stuff) >= 4
                facet = []
                for x in stuff[1:]:
                    facet.append(int(x.split("/")[0]) - 1)
                facets.append(vector_of_unsigned(facet))
    f.close()
    nverts = len(verts)
    # for i in xrange(len(facets)):
    #     for j in xrange(len(facets[i])):
    #         assert facets[i][j] < nverts
    poly = Polyhedron(vector_of_Vector3d(verts), vector_of_vector_of_unsigned(facets))
    return poly

#-------------------------------------------------------------------------------
# Write an OBJ (vertex-facet labeled) shape file from a polyhedron.
# Generalized for polygons now as well.
# NOTE!!  These things seem to count from 1 rather than 0 for indices!
#-------------------------------------------------------------------------------
def writePolyhedronOBJ(poly, filename, forceTriangles=False):
    f = open(filename, "w")
    verts = poly.vertices
    facets = poly.facets
    for v in verts:
        f.write("v %g %g %g\n" % (v.x, v.y, v.z))

    if isinstance(poly, Polygon):
        # There could be multiple loops in a Spheral polygon, so break
        # those up into individual facets for obj output.
        # First, sort topologically so the loops are contiguous
        pairs = [(fac.ipoint1, fac.ipoint2) for fac in facets]
        for i in xrange(len(pairs) - 1):
            v1 = pairs[i][1]
            j = i + 1
            while j < len(pairs) and v1 != pairs[j][0]:
                j += 1
            if j < len(pairs):
                pairs[i+1], pairs[j] = pairs[j], pairs[i+1]
        
        # Now we can write the facets.
        f.write("f")
        i = 0
        for i in xrange(len(pairs)):
            f.write(" %i" % (pairs[i][0] + 1))
            if ((i + 1) < len(pairs) and
                pairs[i][1] != pairs[i+1][0]):
                f.write("\nf")
        f.write("\n")

    else:
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
    verts = []
    for i in xrange(nv):
        line = f.readline()
        stuff = line.split()
        assert len(stuff) == 3
        verts.append(Vector3d(float(stuff[0]), float(stuff[1]), float(stuff[2])))
    
    # Read the face vertex indices.
    facets = []
    for i in xrange(nf):
        line = f.readline()
        stuff = line.split()
        assert len(stuff) >= 4
        n = int(stuff[0])
        facet = []
        for j in xrange(n):
            k = int(stuff[j+1])
            assert k < nv
            facet.append(k)
        facets.append(vector_of_unsigned(facet))
    
    # Build the polyhedron and we're done.
    poly = Polyhedron(vector_of_vector3d(verts), vector_of_vector_of_unsigned(facets))
    return poly


#-------------------------------------------------------------------------------
# Write an OFF shape file from a polyhedron.
# We fudge this a bit and don't write the edges number out.  Hopefully that 
# doesn't bust anything!
#-------------------------------------------------------------------------------
def writePolyhedronOFF(poly, filename):
    f = open(filename, "w")
    verts = poly.vertices
    facets = poly.facets
    f.write("OFF\n")

    # Check for 2d/3d.
    if isinstance(poly, Polygon):
        # There could be multiple loops in a Spheral polygon, so break
        # those up into individual facets for obj output.
        # First, sort topologically so the loops are contiguous
        pairs = [(fac.ipoint1, fac.ipoint2) for fac in facets]
        for i in xrange(len(pairs) - 1):
            v1 = pairs[i][1]
            j = i + 1
            while j < len(pairs) and v1 != pairs[j][0]:
                j += 1
            if j < len(pairs):
                pairs[i+1], pairs[j] = pairs[j], pairs[i+1]
        
        # Now we can extract the faces.
        faces = []
        faces.append([])
        for i in xrange(len(pairs)):
            faces[-1].append(pairs[i][0])
            if (i + 1) < len(pairs) and pairs[i][1] != pairs[i+1][0]:
                faces.append([])

        # Now we can write some header info
        f.write("%i %i %i\n" % (verts.size(), len(faces), 0))

        # Write the vertex coordinates.
        for v in verts:
            f.write("%g %g %g\n" % (v.x, v.y, v.z))

        # Write the faces.
        for face in faces:
            f.write("%i" % len(face))
            for i in face:
                f.write(" %i" % i)
            f.write("\n")

    else:

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

#-------------------------------------------------------------------------------
# Write an STL shape file from a list of polyhedra.
#-------------------------------------------------------------------------------
def writePolyhedraSTL(polys,
                      names,
                      filename):
    assert len(polys) == len(names)
    f = open(filename, "w")
    for name, poly in zip(names, polys):
        verts = poly.vertices
        facets = poly.facets
        f.write("solid %s\n" % name)
        for facet in facets:
            tris = facet.triangles()
            for subfacet in tris:
                ipoints = subfacet.ipoints
                f.write("  facet normal %e %e %e\n" % tuple(subfacet.normal))
                f.write("    outer loop\n")
                for i in ipoints:
                    f.write("      vertex %e %e %e\n" % tuple(verts[i]))
                f.write("    endloop\n")
                f.write("  endfacet\n")
        f.write("endsolid\n")
    f.close()
    return

#-------------------------------------------------------------------------------
# Read an STL polyhedron
#-------------------------------------------------------------------------------
def readPolyhedronSTL(filename):
    mesh = stl.mesh.Mesh.from_file(filename)
    nfacets = len(mesh.v0)
    assert len(mesh.v1) == len(mesh.v2) == nfacets
    points, facets = [], []
    for i in xrange(nfacets):
        points += [Vector3d(*tuple(mesh.v0[i])), Vector3d(*tuple(mesh.v1[i])), Vector3d(*tuple(mesh.v2[i]))]
        facets.append(vector_of_unsigned([3*i, 3*i + 1, 3*i + 2]))
    assert len(points)/3 == nfacets
    assert len(facets) == nfacets
    points = vector_of_Vector3d(points)
    facets = vector_of_vector_of_unsigned(facets)
    poly = Polyhedron(points, facets)
    return poly
