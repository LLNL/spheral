#-------------------------------------------------------------------------------
# JohnsonCookDamageGaussian
#
# Implement the Gaussian distribution for the D1 and D2 fields for the
# Johnson-Cook damage model.
#
# We do this by making a factory function that returns an instance of
# JohnsonCookDamage with the D1 and D2 fields filled in.
#-------------------------------------------------------------------------------
from math import *
import numpy as np
import mpi

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

def JohnsonCookDamageGaussian(nodeList,
                              D1,
                              D2,
                              D3,
                              D4,
                              D5,
                              epsilondot0,
                              Tcrit,
                              sigmamax,
                              efailmin,
                              sigmaD1,
                              sigmaD2,
                              seed = None,
                              domainIndependent = True):
                             
    # What dimension are we?
    assert nodeList.__class__.__name__ in ("SolidNodeList1d", "SolidNodeList2d", "SolidNodeList3d")
    ndim = int(nodeList.__class__.__name__[-2:-1])
    assert ndim in spheralDimensions()

    # Import the approprite bits of Spheral.
    exec("""
from SpheralModules.Spheral.PhysicsSpace import JohnsonCookDamage%(ndim)sd as JohnsonCookDamage
from SpheralModules.Spheral.FieldSpace import ScalarField%(ndim)sd as ScalarField
from SpheralModules.Spheral.DataBaseSpace import DataBase%(ndim)sd as DataBase
from SpheralModules.Spheral import nodeOrdering%(ndim)sd as nodeOrdering
from SpheralModules.Spheral import mortonOrderIndices%(ndim)sd as mortonOrderIndices
""" % {"ndim" : ndim})

    # Prepare the fields for D1 and D2
    fD1 = ScalarField("D1", nodeList, D1)
    fD2 = ScalarField("D2", nodeList, D2)

    # We need to generate the D1 and D2 Fields Gaussian distributed fields.
    # If the a value is 0.0, we take this to mean the corresponding D coefficient is constant.
    if sigmaD1 != 0.0 or sigmaD2 != 0.0:

        # Are we generating domain-independent?
        if domainIndependent:

            # Initialize the random number generator.
            np.random.seed(seed % 2**32)

            # Assign a unique ordering to the nodes so we can step through them
            # in a domain independent manner.
            db = DataBase()
            db.appendNodeList(nodeList)
            keyList = mortonOrderIndices(db);
            orderingList = nodeOrdering(keyList);
            assert orderingList.numFields == 1
            ordering = orderingList[0]
            nglobal = max(0, ordering.max() + 1)

            # Reverse lookup in the ordering.
            order2local = [ordering[i] for i in xrange(nodeList.numInternalNodes)]

            # Generate random numbers in batches.
            iglobal = 0
            nbatch = max(1, nglobal/mpi.procs)
            for iproc in xrange(mpi.procs):
                if iproc == mpi.procs - 1:
                    nlocal = nglobal - iglobal
                else:
                    nlocal = nbatch
                if sigmaD1 != 0.0:
                    D1vals = np.random.normal(D1, sigmaD1, nlocal)
                if sigmaD2 != 0.0:
                    D2vals = np.random.normal(D2, sigmaD2, nlocal)
                for i in xrange(nlocal):
                    try:
                        j = order2local.index(iglobal + i)
                        if sigmaD1 != 0.0:
                            fD1[j] = D1vals[i]
                        if sigmaD2 != 0.0:
                            fD2[j] = D2vals[i]
                    except ValueError:
                        pass
                iglobal += nlocal

        else:

            # In the non-domain independent case we can generate more quickly in parallel.
            procID = mpi.rank
            np.random.seed(((seed + procID)*(seed + procID + 1)/2 + procID) % 2**32)
            if sigmaD1 != 0.0:
                vals = np.random.normal(D1, sigmaD1, nodeList.numInternalNodes)
                for i in xrange(nodeList.numInternalNodes):
                    fD1[i] = vals[i]
            if sigmaD2 != 0.0:
                vals = np.random.normal(D2, sigmaD2, nodeList.numInternalNodes)
                for i in xrange(nodeList.numInternalNodes):
                    fD2[i] = vals[i]

    # Now we build and and return the damage model.
    return JohnsonCookDamage(nodeList,
                             fD1,
                             fD2,
                             D3,
                             D4,
                             D5,
                             epsilondot0,
                             Tcrit,
                             sigmamax,
                             efailmin)
