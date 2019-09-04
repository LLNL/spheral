//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH interpolate.
//------------------------------------------------------------------------------
#include "interpolateCRKSPH.hh"
#include "CRKSPH/CRKSPHUtilities.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "SPH/NodeCoupling.hh"
#include "Utilities/Timer.hh"

extern Timer TIME_interpolateCRKSPH;

namespace Spheral {

using boost::variant;
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace {
//------------------------------------------------------------------------------
// Verify FieldList sizes.
//------------------------------------------------------------------------------
struct FieldListSize: public boost::static_visitor<int> {
  template<typename FieldListType>
  inline
  int operator()(const FieldListType& x) const {
    return x.size();
  }
};

template<typename Dimension> 
inline
bool checkSizes(const vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
                                     FieldList<Dimension, typename Dimension::Vector>,
                                     FieldList<Dimension, typename Dimension::Tensor>,
                                     FieldList<Dimension, typename Dimension::SymTensor>,
                                     FieldList<Dimension, typename Dimension::ThirdRankTensor>>>& fieldLists,
                const size_t numNodeLists) {
  for (const auto& f: fieldLists) {
    if (boost::apply_visitor(FieldListSize(), f) != numNodeLists) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// copyFields
//------------------------------------------------------------------------------
struct CopyFields: public boost::static_visitor<> {
  template<typename FieldListType>
  inline
  void operator()(FieldListType& x) const { x.copyFields(); }
};

//------------------------------------------------------------------------------
// Zero
//------------------------------------------------------------------------------
struct ZeroFields: public boost::static_visitor<> {
  template<typename FieldListType>
  inline
  void operator()(FieldListType& x) const { x.Zero(); }
};

//------------------------------------------------------------------------------
// PrependNameFields
//------------------------------------------------------------------------------
struct PrependNameFields: public boost::static_visitor<> {
  std::string prefix;
  
  PrependNameFields(const std::string prefix_):
    boost::static_visitor<>(),
    prefix(prefix_) {}

  template<typename FieldListType>
  inline
  void operator()(FieldListType& x) const {
    for (auto fieldPtr: x) fieldPtr->name(prefix + fieldPtr->name());
  }
};

//------------------------------------------------------------------------------
// IncrementElement
//------------------------------------------------------------------------------
struct IncrementElement: public boost::static_visitor<> {
  int nodeListi, i, nodeListj, j;
  double mult;

  IncrementElement(int nodeListi_, int i_, int nodeListj_, int j_, double mult_):
    boost::static_visitor<>(),
    nodeListi(nodeListi_),
    i(i_),
    nodeListj(nodeListj_),
    j(j_),
    mult(mult_) {}

  template<typename A, typename B>
  inline
  void operator()(A& x, const B& y) const { VERIFY(false); }

  template<typename FieldListType>
  inline
  void operator()(FieldListType& x, const FieldListType& y) const {
    CHECK2(nodeListi < x.size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    CHECK2(i < x[nodeListi]->size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    CHECK2(nodeListj < y.size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    CHECK2(j < y[nodeListj]->size(),
           "Bad: (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") [("
           << x[nodeListi]->numInternalElements() << " " << x[nodeListi]->size() << ") ("
           << y[nodeListj]->numInternalElements() << " " << y[nodeListj]->size() << ")]");
    x(nodeListi, i) += mult*y(nodeListj, j);
  }
};

}

//------------------------------------------------------------------------------
// interplateCRKSPH
//------------------------------------------------------------------------------
template<typename Dimension>
vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
               FieldList<Dimension, typename Dimension::Vector>,
               FieldList<Dimension, typename Dimension::Tensor>,
               FieldList<Dimension, typename Dimension::SymTensor>,
               FieldList<Dimension, typename Dimension::ThirdRankTensor>>>
interpolateCRKSPH(const vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
                                       FieldList<Dimension, typename Dimension::Vector>,
                                       FieldList<Dimension, typename Dimension::Tensor>,
                                       FieldList<Dimension, typename Dimension::SymTensor>,
                                       FieldList<Dimension, typename Dimension::ThirdRankTensor>>>& fieldLists,
                  const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const FieldList<Dimension, typename Dimension::Scalar>& A,
                  const FieldList<Dimension, typename Dimension::Vector>& B,
                  const FieldList<Dimension, typename Dimension::Tensor>& C,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  const CRKOrder correctionOrder,
                  const TableKernel<Dimension>& W,
                  const NodeCoupling& nodeCoupling) {

  TIME_interpolateCRKSPH.start();

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef vector<variant<FieldList<Dimension, typename Dimension::Scalar>,
                         FieldList<Dimension, typename Dimension::Vector>,
                         FieldList<Dimension, typename Dimension::Tensor>,
                         FieldList<Dimension, typename Dimension::SymTensor>,
                         FieldList<Dimension, typename Dimension::ThirdRankTensor>>> FieldListArray;

  // Pre-conditions.
  const auto numNodeLists = position.size();
  const auto numFieldLists = fieldLists.size();
  REQUIRE(checkSizes<Dimension>(fieldLists, numNodeLists));
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(A.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists or correctionOrder == CRKOrder::ZerothOrder);
  REQUIRE(C.size() == numNodeLists or correctionOrder != CRKOrder::QuadraticOrder);

  // Prepare the result.
  FieldListArray result;
  for (const auto& fieldList: fieldLists) {
    result.push_back(fieldList);
    boost::apply_visitor(CopyFields(), result.back());
    boost::apply_visitor(ZeroFields(), result.back());
    boost::apply_visitor(PrependNameFields("interpolate "), result.back());
  }

  // Walk the FluidNodeLists.
  auto Bi = Vector::zero, Bj = Vector::zero;
  auto Ci = Tensor::zero, Cj = Tensor::zero;
  for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const auto firstGhostNodei = A[nodeListi]->nodeList().firstGhostNode();

    // Iterate over the nodes in this node list.
    for (auto iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const auto i = *iItr;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto& Ai = A(nodeListi, i);
      if (correctionOrder != CRKOrder::ZerothOrder) Bi = B(nodeListi, i);
      if (correctionOrder == CRKOrder::QuadraticOrder) Ci = C(nodeListi, i);

      // Add the self-contribution to each FieldList.
      const auto W0 = W.kernelValue(0.0, Hdeti);
      for (auto k = 0; k < numFieldLists; ++k) {
        boost::apply_visitor(IncrementElement(nodeListi, i, nodeListi, i, weight(nodeListi, i)*W0*Ai), result[k], fieldLists[k]);
      }

      // Neighbors!
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Walk the neighbor nodeLists.
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      
        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const auto firstGhostNodej = A[nodeListj]->nodeList().firstGhostNode();

          // Loop over the neighbors.
          for (auto j: connectivity) {

            // The coupling between these nodes.
            const auto fij = nodeCoupling(nodeListi, i, nodeListj, j);

            // Only proceed if this node pair has not been calculated yet.
            if (fij > 0.0 and connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                                       nodeListj, j,
                                                                       firstGhostNodej)) {

              // Find the effective weights of i->j and j->i.
              // const Scalar wi = fij*2.0*weight(nodeListi, i)*weight(nodeListj, j)/(weight(nodeListi, i) + weight(nodeListj, j));
              // const Scalar wi = fij*0.5*(weight(nodeListi, i) + weight(nodeListj, j));
              // const Scalar wj = wi;
              const auto wi = fij*weight(nodeListi, i);
              const auto wj = fij*weight(nodeListj, j);

              // Get the state for node j.
              const auto& rj = position(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              const auto  Hdetj = Hj.Determinant();
              const auto  Aj = A(nodeListj, j);
              if (correctionOrder != CRKOrder::ZerothOrder) Bj = B(nodeListj, j);
              if (correctionOrder == CRKOrder::QuadraticOrder) Cj = C(nodeListj, j);

              // Node displacement.
              const auto rij = ri - rj;
              const auto etai = Hi*rij;
              const auto etaj = Hj*rij;

              // Kernel weight.
              const auto Wj = CRKSPHKernel(W, correctionOrder,  rij,  etaj, Hdetj, Ai, Bi, Ci);
              const auto Wi = CRKSPHKernel(W, correctionOrder, -rij, -etai, Hdeti, Aj, Bj, Cj);

              // Increment the pair-wise values.
              for (auto k = 0; k < numFieldLists; ++k) {
                boost::apply_visitor(IncrementElement(nodeListi, i, nodeListj, j, wj*Wj), result[k], fieldLists[k]);
                boost::apply_visitor(IncrementElement(nodeListj, j, nodeListi, i, wi*Wi), result[k], fieldLists[k]);
              }

            }
          }
        }
      }
    }
  }

  TIME_interpolateCRKSPH.stop();

  // That's it!
  return result;
}

}
