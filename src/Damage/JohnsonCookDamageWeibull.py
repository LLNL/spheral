#-------------------------------------------------------------------------------
# JohnsonCookDamageWeibull
#
# Implement the Weibull distribution for the D1 and D2 fields for the
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

def JohnsonCookDamageWeibull(nodeList,
                             D1,
                             D2,
                             D3,
                             D4,
                             D5,
                             epsilondot0,
                             Tcrit,
                             sigmamax,
                             efailmin,
                             aD1,
                             bD1,
                             eps0D1,
                             aD2,
                             bD2,
                             eps0D2,
                             seed = None,
                             domainIndependent = True):
                             
    # What dimension are we?
    assert nodeList.__class__.__name__ in ("SolidNodeList1d", "SolidNodeList3d", "SolidNodeList3d")
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
    fD1 = ScalarField("D1", nodeList)
    fD2 = ScalarField("D2", nodeList)

    # We need to generate the D1 and D2 Fields Weibull distributed fields.
    # If the a value is 0.0, we take this to mean the corresponding D coefficient is constant.
    if aD1 != 0.0 or aD2 != 0.0:

        # Are we generating domain-independent?
        if domainIndependent:

            # Initialize the random number generator.
            np.random.seed(seed)

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
                if aD1 != 0.0:
                    D1vals = aD1*np.random.weibull(bD1, nlocal) + eps0D1
                    D2vals = aD2*np.random.weibull(bD2, nlocal) + eps0D2
                    for i in xrange(nlocal):
                        try:
                            j = order2local.index(iglobal + i)
                            if aD1 != 0.0:
                                fD1 = D1vals[i]
                            if aD2 != 0.0:
                                fD2 = D2vals[i]
                        except ValueError:
                            pass
                    iglobal += nlocal

        else:

            # In the non-domain independent case we can generate more quickly in parallel.
            procID = mpi.rank
            np.random.seed((seed + procID)*(seed + procID + 1)/2 + procID)
            if aD1 != 0.0:
                vals = aD1*np.random.weibull(bD1, nodeList.numInternalNodes) + eps0D1
                for i in xrange(nodeList.numInternalNodes):
                    fD1[i] = vals[i]
            if aD2 != 0.0:
                vals = aD2*np.random.weibull(bD2, nodeList.numInternalNodes) + eps0D2
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
