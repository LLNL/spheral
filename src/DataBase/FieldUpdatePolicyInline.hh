#include "Utilities/packElement.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
FieldUpdatePolicy<Dimension, Value>::
FieldUpdatePolicy(std::initializer_list<std::string> depends):
  UpdatePolicyBase<Dimension>(depends) {
}

//------------------------------------------------------------------------------
// Serialize the Field data
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
void
FieldUpdatePolicy<Dimension, Value>::
serializeData(const KeyType& key,
              const State<Dimension>& state,
              std::vector<double>& buf) const {
  const auto& f = state.template field<Value>(key);
  CHECK(f.fixedSizeDataType());
  const std::vector<char> rawbuf = packFieldValues(f);
  const auto ndvals = f.numInternalElements() * f.numValsInDataType();
  const auto nraw = rawbuf.size();
  CHECK(nraw == ndvals*sizeof(double));
  const auto istart = buf.size();
  buf.resize(buf.size() + ndvals);
  CHECK((buf.size() - istart)*sizeof(double) == nraw);
  std::memcpy(&buf[istart], &rawbuf[0], nraw);
  CONTRACT_VAR(ndvals);
}

//------------------------------------------------------------------------------
// Deserialize the Field data
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
size_t
FieldUpdatePolicy<Dimension, Value>::
deserializeData(const KeyType& key,
                const State<Dimension>& state,
                const std::vector<double>& buf,
                const size_t offset) const {
  auto& f = state.template field<Value>(key);
  CHECK(f.fixedSizeDataType());
  const auto ndvals = f.numInternalElements() * f.numValsInDataType();
  CHECK(offset + ndvals <= buf.size());
  const size_t nraw = ndvals * sizeof(double);
  std::vector<char> rawbuf(nraw);
  std::memcpy(&rawbuf[0], &buf[offset], nraw);
  unpackFieldValues(f, rawbuf);
  return offset + ndvals;
}

}
