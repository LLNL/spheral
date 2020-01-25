import mpi
from SpheralCompiledPackages import Vector2d, Vector3d

#-------------------------------------------------------------------------------
# Reject in a polygonal surface
#-------------------------------------------------------------------------------
class PolygonalSurfaceRejecter:

    def __init__(self, surface,
                 interior = True):
        self.surface = surface
        self.interior = interior # Reject interior to surface?
        return

    # New style rejecter
    def accept(self, x, y):
        return self.interior ^ self.surface.contains(Vector2d(x,y))

    # Old style rejecter
    def __call__(self, x, y, m, H):
        n = len(x)
        assert len(y) == n
        assert len(m) == n
        assert len(H) == n

        # We'll take advantage of any available parallelism to split
        # up the containment testing.  The following algorithm is borrowed 
        # from NodeGeneratorBase to divvy up the ID range.
        ndomain0 = n/mpi.procs
        remainder = n % mpi.procs
        assert remainder < mpi.procs
        ndomain = ndomain0
        if mpi.rank < remainder:
            ndomain += 1
        imin = mpi.rank*ndomain0 + min(mpi.rank, remainder)
        imax = imin + ndomain

        # Check our local range of IDs.
        xloc, yloc, zloc, mloc, Hloc = [], [], [], [], []
        localIndices = [i for i in xrange(imin, imax) 
                        if self.accept(x[i], y[i])]

        # Now cull to the interior values.
        xnew, ynew, mnew, Hnew = [], [], [], []
        for iproc in xrange(mpi.procs):
            otherIndices = mpi.bcast(localIndices, iproc)
            for i in otherIndices:
                xnew.append(x[i])
                ynew.append(y[i])
                mnew.append(m[i])
                Hnew.append(H[i])

        # That's it.
        return xnew, ynew, mnew, Hnew

#-------------------------------------------------------------------------------
# Reject in a polyhedral surface
#-------------------------------------------------------------------------------
class PolyhedralSurfaceRejecter:

    def __init__(self, surface,
                 interior = True):
        self.surface = surface
        self.interior = interior # Reject interior to surface?
        return

    # New style rejecter
    def accept(self, x, y, z):
        return self.interior ^ self.surface.contains(Vector3d(x,y,z))

    # Old style rejecter
    def __call__(self, x, y, z, m, H):
        n = len(x)
        assert len(y) == n
        assert len(z) == n
        assert len(m) == n
        assert len(H) == n

        # We'll take advantage of any available parallelism to split
        # up the containment testing.  The following algorithm is borrowed 
        # from NodeGeneratorBase to divvy up the ID range.
        ndomain0 = n/mpi.procs
        remainder = n % mpi.procs
        assert remainder < mpi.procs
        ndomain = ndomain0
        if mpi.rank < remainder:
            ndomain += 1
        imin = mpi.rank*ndomain0 + min(mpi.rank, remainder)
        imax = imin + ndomain

        # Check our local range of IDs.
        xloc, yloc, zloc, mloc, Hloc = [], [], [], [], []
        localIndices = [i for i in xrange(imin, imax) 
                        if self.accept(x[i], y[i], z[i])]

        # Now cull to the interior values.
        xnew, ynew, znew, mnew, Hnew = [], [], [], [], []
        for iproc in xrange(mpi.procs):
            otherIndices = mpi.bcast(localIndices, iproc)
            for i in otherIndices:
                xnew.append(x[i])
                ynew.append(y[i])
                znew.append(z[i])
                mnew.append(m[i])
                Hnew.append(H[i])

        # That's it.
        return xnew, ynew, znew, mnew, Hnew
