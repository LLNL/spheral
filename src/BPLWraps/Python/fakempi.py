#-------------------------------------------------------------------------------
# A fake version of the pyMPI mpi module for use in serial runs.
#-------------------------------------------------------------------------------
class fakempi:

    def __init__(self):
        print "Invoking fake mpi module."
        self.rank = 0
        self.procs = 1
        self.MIN = -1
        self.MAX = -2
        self.SUM = -3
        return

    def reduce(self, var, op):
        return var

    def allreduce(self, var, op):
        return var

    def gather(self, obj, op):
        return obj

    def bcast(self, obj, root=0):
        return obj

    def barrier(self):
        return
