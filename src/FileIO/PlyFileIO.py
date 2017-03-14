from Spheral3d import *

def PolyFromPly(filename):

    points = 0
    faces  = 0

    f = open(filename,'r')
    for line in f:
        data = line.split()
        if data[0] == "element":
            if data[1] == "vertex":
                points = int(data[2])
            if data[1] == "face":
                faces = int(data[2])
    f.close()

    print "Got %d points in %d faces" % (points,faces)

    vertices = []
    indices  = []

    start = 100
    i = 0
    f = open(filename,'r')
    for line in f:
        data = line.split()
        if i >= start:
            if len(vertices) == points:
                indices.append(data[1:])
            else:
                vertices.append(data[:3])
        elif data[0] == "end_header":
            #print "found end_header at line %d, starting read at line %d" % (i,i+1)
            start = i + 1
        i += 1
    f.close()
    for i in xrange(len(vertices)):
        for j in xrange(len(vertices[i])):
            vertices[i][j] = float(vertices[i][j])
    for i in xrange(len(indices)):
        for j in xrange(len(indices[i])):
            indices[i][j] = int(indices[i][j])

    vvert = vector_of_Vector()
    vindx = vector_of_vector_of_unsigned()

    for i in xrange(len(vertices)):
        vvert.append(Vector(vertices[i]))
    for i in xrange(len(indices)):
        vuns = vector_of_unsigned()
        for j in xrange(len(indices[i])):
            vuns.append(indices[i][j])
        vindx.append(vuns)

    shape = Polyhedron(vvert,vindx)
    return shape

