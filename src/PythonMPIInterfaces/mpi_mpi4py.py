#-------------------------------------------------------------------------------
# mpi
#
# This module reproduces the pyMPI interface using mpi4py.
#-------------------------------------------------------------------------------
import sys
from SpheralTestUtilities import globalFrame
from mpi4py import MPI

#-------------------------------------------------------------------------------
# Communicator and geometry.
#-------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
procs = comm.Get_size()

#-------------------------------------------------------------------------------
# Define the operations.
#-------------------------------------------------------------------------------
MIN = MPI.MIN
MAX = MPI.MAX
SUM = MPI.SUM

#-------------------------------------------------------------------------------
# Prepare files to keep the stdout and stderr streams in.
# The pyMPI defaults are only rank 0 writes stdout, but all
# processes write stderr.
#-------------------------------------------------------------------------------
globalscope = globalFrame().f_globals
if rank > 0:
    exec("""
import sys
__mpi_stdoutfile__ = open("/dev/null", "w")
sys.stdout = __mpi_stdoutfile__
""", globalscope)

#-------------------------------------------------------------------------------
# send
#-------------------------------------------------------------------------------
def send(obj, dest=0, tag=100):
    comm.send(obj=obj, dest=dest, tag=tag)

#-------------------------------------------------------------------------------
# recv
#-------------------------------------------------------------------------------
def recv(source=0, tag=100, obj=None):
    return (comm.recv(obj=obj, source=source, tag=tag), )

#-------------------------------------------------------------------------------
# isend
#-------------------------------------------------------------------------------
def isend(obj, dest=0, tag=100):
    comm.isend(obj=obj, dest=dest, tag=tag)

#-------------------------------------------------------------------------------
# reduce
#-------------------------------------------------------------------------------
def reduce(obj, op=SUM, root=0):
    return comm.reduce(obj, op=op, root=root)

#-------------------------------------------------------------------------------
# allreduce
#-------------------------------------------------------------------------------
def allreduce(obj, op=SUM):
    return comm.allreduce(obj, op=op)

#-------------------------------------------------------------------------------
# gather
#-------------------------------------------------------------------------------
def gather(obj, root=0):
    return comm.gather(obj, root=root)

#-------------------------------------------------------------------------------
# allgather
#-------------------------------------------------------------------------------
def allgather(obj):
    return comm.allgather(obj)

#-------------------------------------------------------------------------------
# bcast
#-------------------------------------------------------------------------------
def bcast(obj, root=0):
    return comm.bcast(obj, root=root)

#-------------------------------------------------------------------------------
# barrier
#-------------------------------------------------------------------------------
def barrier():
    comm.barrier()

#-------------------------------------------------------------------------------
# synchronizeQueuedOutput
#-------------------------------------------------------------------------------
def synchronizeQueuedOutput(stdoutfile = None,
                            stderrfile = None):
    if stdoutfile == None:
        exec("sys.stdout = sys.__stdout__", globalscope)
    else:
        exec("__mpi_stdoutfile__ = open(%s, 'w'); sys.stdout = __mpi_stdoutfile__" % stdoutfile,
             globalscope)

    if stderrfile == None:
        exec("sys.stderr = sys.__stderr__", globalscope)
    else:
        exec("__mpi_stderrfile__ = open(%s, 'w'); sys.stderr = __mpi_stderrfile__" % stderrfile,
             globalscope)

    return
