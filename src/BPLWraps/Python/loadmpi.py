#-------------------------------------------------------------------------------
# Helper function to safely get the mpi module, returning a dummy when we're
# using non-mpi serial python
#-------------------------------------------------------------------------------
import fakempi
def loadmpi():
    try:
        import mpi
        procID = mpi.rank
        nProcs = mpi.procs
    except:
        mpi = fakempi.fakempi()
        procID = 0
        nProcs = 1

    return mpi, procID, nProcs

