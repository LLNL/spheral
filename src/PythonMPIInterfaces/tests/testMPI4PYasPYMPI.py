#ATS:test(SELF, np=2, label="Test compatibility interface faking pyMPI using mpi4py (2 proc)")
#ATS:test(SELF, np=4, label="Test compatibility interface faking pyMPI using mpi4py (4 proc)")

from math import *
import unittest
from SpheralTestUtilities import *

import mpi
from mpi4py import MPI

#-------------------------------------------------------------------------------
# Test our interface class reproducing the pyMPI interface using mpi4py.
#-------------------------------------------------------------------------------
class MPI4PYasPYMPI(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.comm = MPI.COMM_WORLD
        return

    #---------------------------------------------------------------------------
    # tearDown
    #---------------------------------------------------------------------------
    def tearDown(self):
        return

    #---------------------------------------------------------------------------
    # check -- utility to stop all process if there's a failure
    #---------------------------------------------------------------------------
    def check(self, cond):
        x = 0
        if cond:
            x = 1
        globalCheck = self.comm.allreduce(x, op=MPI.MIN)
        assert globalCheck == 1

    #---------------------------------------------------------------------------
    # rank
    #---------------------------------------------------------------------------
    def testRank(self):
        self.check(mpi.rank == self.comm.Get_rank())

    #---------------------------------------------------------------------------
    # procs
    #---------------------------------------------------------------------------
    def testProcs(self):
        self.check(mpi.procs == self.comm.Get_size())

    #---------------------------------------------------------------------------
    # generic reduce test
    #---------------------------------------------------------------------------
    def reduceObjImpl(self, x, op, ans):
        for root in xrange(mpi.procs):
            globalx = mpi.reduce(x, op=op, root=root)
            ok = True
            if mpi.rank == root:
                ok = (globalx == ans)
            else:
                ok = (globalx is None)
            if not ok:
                print globalx, ans
            self.check(ok)

    #---------------------------------------------------------------------------
    # generic allreduce test
    #---------------------------------------------------------------------------
    def allreduceObjImpl(self, x, op, ans):
        globalx = mpi.allreduce(x, op=op)
        ok = (globalx == ans)
        if not ok:
            sys.stderr.write("failure:  %s != %s\n" % (str(globalx), str(ans)))
        self.check(ok)

    #---------------------------------------------------------------------------
    # reduce float MIN
    #---------------------------------------------------------------------------
    def testReduceFloatMin(self):
        x = 10.0*mpi.rank
        self.reduceObjImpl(x, mpi.MIN, 0.0)

    #---------------------------------------------------------------------------
    # reduce float MAX
    #---------------------------------------------------------------------------
    def testReduceFloatMax(self):
        x = 10.0*mpi.rank
        self.reduceObjImpl(x, mpi.MAX, 10.0*(mpi.procs - 1))

    #---------------------------------------------------------------------------
    # reduce float SUM
    #---------------------------------------------------------------------------
    def testReduceFloatSum(self):
        x = 10.0*mpi.rank
        self.reduceObjImpl(x, mpi.SUM, 10.0*sum(range(mpi.procs)))

    #---------------------------------------------------------------------------
    # reduce list SUM
    #---------------------------------------------------------------------------
    def testReduceListSum(self):
        x = range(10*mpi.rank, 10*mpi.rank + 10)
        self.reduceObjImpl(x, mpi.SUM, range(10*mpi.procs))

    #---------------------------------------------------------------------------
    # allreduce float MIN
    #---------------------------------------------------------------------------
    def testAllreduceFloatMin(self):
        x = 10.0*mpi.rank
        self.allreduceObjImpl(x, mpi.MIN, 0.0)

    #---------------------------------------------------------------------------
    # allreduce float MAX
    #---------------------------------------------------------------------------
    def testAllreduceFloatMax(self):
        x = 10.0*mpi.rank
        self.allreduceObjImpl(x, mpi.MAX, 10.0*(mpi.procs - 1))

    #---------------------------------------------------------------------------
    # allreduce float SUM
    #---------------------------------------------------------------------------
    def testAllreduceFloatSum(self):
        x = 10.0*mpi.rank
        self.allreduceObjImpl(x, mpi.SUM, 10.0*sum(range(mpi.procs)))

    #---------------------------------------------------------------------------
    # gather int
    #---------------------------------------------------------------------------
    def testGatherInt(self):
        x = 10*mpi.rank
        ans = range(0, 10*mpi.procs, 10)
        for root in xrange(mpi.procs):
            globalx = mpi.gather(x, root = root)
            if mpi.rank == root:
                self.check(globalx == ans)
            else:
                self.check(globalx == None)

    #---------------------------------------------------------------------------
    # allgather int
    #---------------------------------------------------------------------------
    def testAllgatherInt(self):
        x = 10*mpi.rank
        ans = range(0, 10*mpi.procs, 10)
        globalx = mpi.allgather(x)
        self.check(globalx == ans)

    #---------------------------------------------------------------------------
    # allgather list(int)
    #---------------------------------------------------------------------------
    def testAllgatherInt(self):
        x = range(10*mpi.rank, 10*mpi.rank + 10)
        ans = []
        for i in xrange(mpi.procs):
            ans.append(range(10*i, 10*i + 10))
        globalx = mpi.allgather(x)
        self.check(globalx == ans)

    #---------------------------------------------------------------------------
    # blocking send/recv
    #---------------------------------------------------------------------------
    def testBlockSend(self):
        for sendProc in xrange(mpi.procs):
            if mpi.rank == sendProc:
                for j in xrange(mpi.procs):
                    if j != mpi.rank:
                        obj = 10*mpi.rank + 1
                        mpi.send(obj, dest=j, tag=100)
            else:
                obj = mpi.recv(sendProc, 100)[0]
                assert obj == 10*sendProc + 1
                
    #---------------------------------------------------------------------------
    # non-blocking send/blocking recv
    #---------------------------------------------------------------------------
    def testNonBlockSend(self):
        for sendProc in xrange(mpi.procs):
            if mpi.rank == sendProc:
                for j in xrange(mpi.procs):
                    if j != mpi.rank:
                        obj = 10*mpi.rank + 1
                        mpi.isend(obj, dest=j, tag=100)
            else:
                obj = mpi.recv(sendProc, 100)[0]
                assert obj == 10*sendProc + 1
                
#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
