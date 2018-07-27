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
import np
import mpi

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

def JohnsonCookDamageWeibull(nodeList,
                             D1,
                             D2,
                             D3,
                             D4,
                             D5,
                             aD1,
                             bD1,
                             eps0D1,
                             aD2,
                             bD2,
                             eps0D2,
                             epsilondot0,
                             Tcrit,
                             sigmamax,
                             efailmin,
                             seed,
                             domainIndependent):
                             
    # What dimension are we?
    assert nodeList.__name__ in ("SolidNodeList1d", "SolidNodeList3d", "SolidNodeList3d")
    ndim = int(nodeList.__name__[-3:-2])
    assert ndim in spheralDimensions()

    # Import the approprite bits of Spheral.
    exec("""
from SpheralModules.Spheral.PhysicsSpace import JohnsonCookDamage%(ndim)sd as JohnsonCookDamage
from SpheralModules.Spheral.FieldSpace import ScalarField%(ndim)sd as ScalarField
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

      // Construct a random number generator.
      // C++11 provides a Weibull distribution natively.
      std::mt19937 gen(seed);
      std::weibull_distribution<> d1(aD1, bD1), d2(aD2, bD2);

      // First, assign a unique ordering to the nodes so we can step through them
      // in a domain independent manner.
      typedef KeyTraits::Key Key;
      DataBase<Dimension> db;
      db.appendNodeList(nodeList);
      auto keyList = mortonOrderIndices(db);
      auto orderingList = nodeOrdering(keyList);
      CHECK(orderingList.numFields() == 1);
      auto& ordering = *orderingList[0];
      const auto n = std::max(0, ordering.max());  // Note this is the global number of nodes.

      // Reverse lookup in the ordering.
      unordered_map<unsigned, unsigned> order2local;
      for (auto i = 0; i != nodeList.numInternalNodes(); ++i) order2local[ordering[i]] = i;

      // Walk the global number of node in Morton order.
      for (auto iorder = 0; iorder < n + 1; ++iorder) {

        // Is this one of our nodes?
        auto itr = order2local.find(iorder);
        if (itr != order2local.end()) {

          // This node is on this domain, so assign D1 and D2.
          const auto i = itr->second;
          CHECK(i < nodeList.numInternalNodes());
          if (aD1 != 0.0) mD1[i] = D1*(d1(gen) + eps0D1);
          if (aD2 != 0.0) mD2[i] = D2*(d2(gen) + eps0D2);

        } else {

          // Otherwise we just spin the random generators to keep all domains in sync.
          d1(gen);
          d2(gen);

        }
      }

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
