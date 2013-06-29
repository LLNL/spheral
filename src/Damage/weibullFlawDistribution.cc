//------------------------------------------------------------------------------
// Construct a flaw distribution for a set of nodes according to a Weibull 
// (power-law) distribution.
//------------------------------------------------------------------------------
#include <set>
#include <algorithm>
#include "boost/unordered_map.hpp"

#include "weibullFlawDistribution.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/mortonOrderIndicies.hh"
#include "Utilities/nodeOrdering.hh"
#include "NodeList/FluidNodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/allReduce.hh"

#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>

namespace Spheral {
namespace PhysicsSpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::FluidNodeList;
using DataBaseSpace::DataBase;
using boost::unordered_map;

//------------------------------------------------------------------------------
// This version uses the Benz-Asphaug algorithm, stepping up deterministically
// in flaw energy based on a minimum chosen from the simulation volume.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, vector<double> >
weibullFlawDistributionBenzAsphaug(double volume,
                                   const double volumeStretchFactor,
                                   const unsigned seed,
                                   const double kWeibull,
                                   const double mWeibull,
                                   const FluidNodeList<Dimension>& nodeList,
                                   const int minFlawsPerNode,
                                   const int minTotalFlaws) {

  // Pre-conditions.
  REQUIRE(volume >= 0.0);
  REQUIRE(volumeStretchFactor >= 1.0);
  REQUIRE(kWeibull >= 0.0);
  REQUIRE(mWeibull > 0.0);
  REQUIRE(minFlawsPerNode > 0);
  REQUIRE(minTotalFlaws > 0);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Prepare the result.
  Field<Dimension, vector<double> > flaws("Weibull flaw distribution",
                                          nodeList);

  // Construct unique global IDs for all nodes in the NodeList.
  const int n = max(1, numGlobalNodes(nodeList));
  const Field<Dimension, int> globalIDs = globalNodeIDs(nodeList);

  // Prepare a table to faciliate looking local IDs from global.
  unordered_map<unsigned, unsigned> global2local;
  for (unsigned i = 0; i != nodeList.numInternalNodes(); ++i) global2local[globalIDs(i)] = i;
  CHECK(global2local.size() == nodeList.numInternalNodes());

  // Prepare an int per *each* node, so that each process can keep track of how many
  // flaws are seeded globally (avoiding communication at the expense of memory).
  vector<int> numFlawsPerNode((size_t) n, 0);

  // Identify the rank and number of domains.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // Only proceed if there are nodes to initialize!
  if (n > 0) {

    // If the user did not speicify a volume, we compute it from the information
    // in the NodeList.
    if (volume == 0.0) {
      const Field<Dimension, Scalar>& mass = nodeList.mass();
      const Field<Dimension, Scalar>& rho = nodeList.massDensity();
      for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
        CHECK(rho(i) > 0.0);
        volume += mass(i)/rho(i);
      }
#ifdef USE_MPI
      volume = allReduce(volume, MPI_SUM, Communicator::communicator());
#endif
    }
    volume = std::max(volume, 1e-100);
    CHECK(volume > 0.0);

    // Compute the minimum (starting) failure strain.
    const double mInv = 1.0/(mWeibull + 1.0e-50);
    const double epsMin = pow(kWeibull*volume*volumeStretchFactor, -mInv);
    CHECK(epsMin > 0.0);

    // Construct a random number generator.
    typedef boost::mt19937 base_generator_type;
    base_generator_type basegen(seed);
    boost::uniform_01<base_generator_type> generator(basegen);

    // Loop and initialize flaws until:
    // a) every node has the minimum number of flaws per node, and
    // b) we meet the minimum number of total flaws.
    int numCompletedNodes = 0;
    int ienergy = 1;
    while ((numCompletedNodes < n) || (ienergy <= minTotalFlaws)) {

      // Randomly select a global node.
      const int iglobal = int(generator() * n);
      CHECK(iglobal >= 0 && iglobal < n);

      // Increment the number of flaws for this node, and check if this
      // completes this node.
      ++numFlawsPerNode[iglobal];
      if (numFlawsPerNode[iglobal] == minFlawsPerNode) ++numCompletedNodes;

      // Is this node one of ours?
      const typename unordered_map<unsigned, unsigned>::const_iterator itr = global2local.find(iglobal);
      if (itr != global2local.end()) {

        const unsigned i = itr->second;
        CHECK(i < nodeList.numInternalNodes());

        // The activation energy.
        const double epsij = epsMin * pow(ienergy*volumeStretchFactor, mInv);

        // Add a flaw with this activation energy to this node.
        flaws(i).push_back(epsij);
      }

      // Increment the energy multiplier.
      ++ienergy;
    }

    // Sort the flaws on each node by energy.
    size_t minNumFlaws = INT_MAX;
    size_t maxNumFlaws = 0;
    size_t totalNumFlaws = 0;
    double epsMax = 0.0;
    double sumFlaws = 0.0;
    for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
      sort(flaws(i).begin(), flaws(i).end());
      minNumFlaws = min(minNumFlaws, flaws(i).size());
      maxNumFlaws = max(maxNumFlaws, flaws(i).size());
      totalNumFlaws += flaws(i).size();
      epsMax = max(epsMax, flaws(i).back());
      for (int j = 0; j != flaws(i).size(); ++j) sumFlaws += flaws(i)[j];
    }

    // Prepare some diagnostic output.
    {
#ifdef USE_MPI
      {
        size_t minNumFlaws0 = minNumFlaws;
        size_t maxNumFlaws0 = maxNumFlaws;
        size_t totalNumFlaws0 = totalNumFlaws;
        double epsMax0 = epsMax;
        double sumFlaws0 = sumFlaws;
        MPI_Allreduce(&minNumFlaws0, &minNumFlaws, 1, MPI_UNSIGNED, MPI_MIN, Communicator::communicator());
        MPI_Allreduce(&maxNumFlaws0, &maxNumFlaws, 1, MPI_UNSIGNED, MPI_MAX, Communicator::communicator());
        MPI_Allreduce(&totalNumFlaws0, &totalNumFlaws, 1, MPI_UNSIGNED, MPI_SUM, Communicator::communicator());
        MPI_Allreduce(&epsMax0, &epsMax, 1, MPI_DOUBLE, MPI_MAX, Communicator::communicator());
        MPI_Allreduce(&sumFlaws0, &sumFlaws, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
      }
#endif
      if (procID == 0) {
        cerr << "weibullFlawDistributionBenzAsphaug: Min num flaws per node: " << minNumFlaws << endl
             << "                                    Max num flaws per node: " << maxNumFlaws << endl
             << "                                    Total num flaws       : " << totalNumFlaws << endl
             << "                                    Avg flaws per node    : " << totalNumFlaws / n << endl
             << "                                    Min flaw strain       : " << epsMin << endl
             << "                                    Max flaw strain       : " << epsMax << endl
             << "                                    Avg node failure      : " << sumFlaws / n << endl;
      }
    }

  }

  // That's it.
  BEGIN_CONTRACT_SCOPE;
  {
    for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
      ENSURE(flaws(i).size() >= minFlawsPerNode);
      for (vector<double>::const_iterator itr = flaws(i).begin() + 1;
           itr != flaws(i).end();
           ++itr) ENSURE(*itr >= *(itr - 1));
    }
  }
  END_CONTRACT_SCOPE;

  return flaws;
}

//------------------------------------------------------------------------------
// This version uses my own algorithm, stochastically seeding flaws in the range
// [0, epsmax] where epsmax is chosen per node based on the nodal volume.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, vector<double> >
weibullFlawDistributionOwen(const unsigned seed,
                            const double kWeibull,
                            const double mWeibull,
                            const FluidNodeList<Dimension>& nodeList,
                            const int numFlawsPerNode,
                            const double volumeMultiplier) {

  // Pre-conditions.
  REQUIRE(kWeibull >= 0.0);
  REQUIRE(mWeibull > 0.0);
  REQUIRE(numFlawsPerNode > 0);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef KeyTraits::Key Key;

  // Prepare the result.
  Field<Dimension, vector<double> > flaws("Weibull flaw distribution",
                                          nodeList);

  // Assign a unique ordering to the nodes so we can step through them
  // in a domain independent manner.
  DataBase<Dimension> db;
  db.appendNodeList(const_cast<FluidNodeList<Dimension>&>(nodeList));
  FieldList<Dimension, Key> keyList = mortonOrderIndicies(db);
  FieldList<Dimension, int> orderingList = nodeOrdering(keyList);
  CHECK(orderingList.numFields() == 1);
  Field<Dimension, int>& ordering = *orderingList[0];
  const int n = std::max(0, ordering.max());

  // Reverse lookup in the ordering.
  unordered_map<unsigned, unsigned> order2local;
  for (unsigned i = 0; i != nodeList.numInternalNodes(); ++i) order2local[ordering[i]] = i;

  // Identify the rank and number of domains.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // State for this NodeList.
  const Field<Dimension, Scalar>& mass = nodeList.mass();
  const Field<Dimension, Scalar>& rho = nodeList.massDensity();

  // Construct a random number generator.
  typedef boost::mt19937 base_generator_type;
  base_generator_type basegen(seed);
  boost::uniform_01<base_generator_type> generator(basegen);

  // Iterate over the nodes.
  const double mInv = 1.0/mWeibull;
  for (int iorder = 0; iorder != n + 1; ++iorder) {

    // Is this one of our nodes?
    typename unordered_map<unsigned, unsigned>::const_iterator itr = order2local.find(iorder);
    if (itr != order2local.end()) {

      // We have the node!
      const unsigned i = itr->second;
      CHECK(i < nodeList.numInternalNodes());
      CHECK(rho(i) > 0.0);
      const double Vi = mass(i)/rho(i) * volumeMultiplier;
      CHECK(Vi > 0.0);
      const double epsmaxi = pow(kWeibull*Vi, -mInv);
      CHECK(epsmaxi > 0.0);
      const double Bi = pow(epsmaxi, mWeibull);
      CHECK(Bi > 0.0);

      // Seed flaws on the node.
      for (int j = 0; j != numFlawsPerNode; ++j) {
        flaws(i).push_back(pow(Bi * generator(), mInv));
        CHECK(flaws(i).back() >= 0.0 and flaws(i).back() <= epsmaxi);
      }
      CHECK(flaws(i).size() == numFlawsPerNode);

    } else {
      // Other domains just cycle the random number generator so that
      // we can be domain decomposition independent.
      for (int j = 0; j != numFlawsPerNode; ++j) {
        double tmp = generator();
      }
    }
  }

  // Sort the flaws on each node by energy.
  size_t minNumFlaws = std::numeric_limits<int>::max();
  size_t maxNumFlaws = 0;
  size_t totalNumFlaws = 0;
  double epsMin = std::numeric_limits<double>::max();
  double epsMax = std::numeric_limits<double>::min();
  double sumFlaws = 0.0;
  for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
    sort(flaws(i).begin(), flaws(i).end());
    minNumFlaws = min(minNumFlaws, flaws(i).size());
    maxNumFlaws = max(maxNumFlaws, flaws(i).size());
    totalNumFlaws += flaws(i).size();
    epsMin = min(epsMin, flaws(i).front());
    epsMax = max(epsMax, flaws(i).back());
    for (int j = 0; j != flaws(i).size(); ++j) sumFlaws += flaws(i)[j];
  }

  // Prepare some diagnostic output.
  if (n > 0) {
#ifdef USE_MPI
    {
      size_t minNumFlaws0 = minNumFlaws;
      size_t maxNumFlaws0 = maxNumFlaws;
      size_t totalNumFlaws0 = totalNumFlaws;
      double epsMin0 = epsMin;
      double epsMax0 = epsMax;
      double sumFlaws0 = sumFlaws;
      MPI_Allreduce(&minNumFlaws0, &minNumFlaws, 1, MPI_UNSIGNED, MPI_MIN, Communicator::communicator());
      MPI_Allreduce(&maxNumFlaws0, &maxNumFlaws, 1, MPI_UNSIGNED, MPI_MAX, Communicator::communicator());
      MPI_Allreduce(&totalNumFlaws0, &totalNumFlaws, 1, MPI_UNSIGNED, MPI_SUM, Communicator::communicator());
      MPI_Allreduce(&epsMin0, &epsMin, 1, MPI_DOUBLE, MPI_MIN, Communicator::communicator());
      MPI_Allreduce(&epsMax0, &epsMax, 1, MPI_DOUBLE, MPI_MAX, Communicator::communicator());
      MPI_Allreduce(&sumFlaws0, &sumFlaws, 1, MPI_DOUBLE, MPI_SUM, Communicator::communicator());
    }
#endif
    if (procID == 0) {
      cerr << "weibullFlawDistributionOwen: Min num flaws per node: " << minNumFlaws << endl
           << "                             Max num flaws per node: " << maxNumFlaws << endl
           << "                             Total num flaws       : " << totalNumFlaws << endl
           << "                             Avg flaws per node    : " << totalNumFlaws / n << endl
           << "                             Min flaw strain       : " << epsMin << endl
           << "                             Max flaw strain       : " << epsMax << endl
           << "                             Avg node failure      : " << sumFlaws / n << endl;
    }
  }

  // That's it.
  BEGIN_CONTRACT_SCOPE;
  {
    for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
      ENSURE(flaws(i).size() == numFlawsPerNode);
      for (vector<double>::const_iterator itr = flaws(i).begin() + 1;
           itr != flaws(i).end();
           ++itr) ENSURE(*itr >= *(itr - 1));
    }
  }
  END_CONTRACT_SCOPE;

  return flaws;
}

}
}

