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

//------------------------------------------------------------------------------
// IncrementFieldList
//------------------------------------------------------------------------------
struct IncrementFieldList: public boost::static_visitor<> {

  template<typename A, typename B>
  inline
  void operator()(A& x, const B& y) const { VERIFY(false); }

  template<typename FieldListType>
  inline
  void operator()(FieldListType& x, const FieldListType& y) const {
    CHECK(x.size() == y.size());
    x += y;
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
                  const RKOrder correctionOrder,
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
  REQUIRE(B.size() == numNodeLists or correctionOrder == RKOrder::ZerothOrder);
  REQUIRE(C.size() == numNodeLists or correctionOrder != RKOrder::QuadraticOrder);

  // Prepare the result.
  FieldListArray result;
  for (const auto& fieldList: fieldLists) {
    result.push_back(fieldList);
    boost::apply_visitor(CopyFields(), result.back());
    boost::apply_visitor(ZeroFields(), result.back());
    boost::apply_visitor(PrependNameFields("interpolate "), result.back());
  }

  // Walk all the interacting pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {

    // Thread private result
    FieldListArray localResult;
    #pragma omp critical
    for (const auto& fieldList: fieldLists) {
      localResult.push_back(fieldList);
      boost::apply_visitor(CopyFields(), localResult.back());
      boost::apply_visitor(ZeroFields(), localResult.back());
      boost::apply_visitor(PrependNameFields("local interpolate "), localResult.back());
    }

    auto Bi = Vector::zero, Bj = Vector::zero;
    auto Ci = Tensor::zero, Cj = Tensor::zero;
    int i, j, nodeListi, nodeListj;

#pragma omp for
    for (auto kk = 0; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto& Ai = A(nodeListi, i);
      if (correctionOrder != RKOrder::ZerothOrder) Bi = B(nodeListi, i);
      if (correctionOrder == RKOrder::QuadraticOrder) Ci = C(nodeListi, i);

      // The coupling between these nodes.
      const auto fij = nodeCoupling(nodeListi, i, nodeListj, j);
      if (fij > 0.0) {

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
        if (correctionOrder != RKOrder::ZerothOrder) Bj = B(nodeListj, j);
        if (correctionOrder == RKOrder::QuadraticOrder) Cj = C(nodeListj, j);

        // Node displacement.
        const auto rij = ri - rj;
        const auto etai = Hi*rij;
        const auto etaj = Hj*rij;

        // Kernel weight.
        const auto Wj = CRKSPHKernel(W, correctionOrder,  rij,  etaj, Hdetj, Ai, Bi, Ci);
        const auto Wi = CRKSPHKernel(W, correctionOrder, -rij, -etai, Hdeti, Aj, Bj, Cj);

        // Increment the pair-wise values.
        for (auto k = 0; k < numFieldLists; ++k) {
          boost::apply_visitor(IncrementElement(nodeListi, i, nodeListj, j, wj*Wj), localResult[k], fieldLists[k]);
          boost::apply_visitor(IncrementElement(nodeListj, j, nodeListi, i, wi*Wi), localResult[k], fieldLists[k]);
        }
      }
    }

    // Merge the local to global result
#pragma omp critical
    {
      for (auto k = 0; k < numFieldLists; ++k) {
        boost::apply_visitor(IncrementFieldList(), result[k], localResult[k]);
      }
    }
  }

  // Add the self contribution.
  auto Bi = Vector::zero;
  auto Ci = Tensor::zero;
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = A[nodeListi]->nodeList().numInternalNodes();
#pragma omp parallel for
    for (auto i = 0; i < n; ++i) {

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto& Ai = A(nodeListi, i);
      if (correctionOrder != RKOrder::ZerothOrder) Bi = B(nodeListi, i);
      if (correctionOrder == RKOrder::QuadraticOrder) Ci = C(nodeListi, i);

      // Add the self-contribution to each FieldList.
      const auto W0 = W.kernelValue(0.0, Hdeti);
      for (auto k = 0; k < numFieldLists; ++k) {
        boost::apply_visitor(IncrementElement(nodeListi, i, nodeListi, i, weight(nodeListi, i)*W0*Ai), result[k], fieldLists[k]);
      }
    }
  }

  TIME_interpolateCRKSPH.stop();

  // That's it!
  return result;
}

}
