#ifndef TENSORTYPE_DEFINED
#define TENSORTYPE_DEFINED

namespace Spheral {

struct FullTensorType;
struct SymmetricTensorType;

struct SymmetricTensorType {
  typedef FullTensorType OtherTensorType;
};
struct FullTensorType {
  typedef SymmetricTensorType OtherTensorType;
};

// enum TensorType {
//   FullTensor = 0,
//   SymmetricTensor = 1
// };

// // The following bit of chicanery is designed to let a (Sym,Full)Tensor
// // know about the other (Full,Sym)Tensor type and grant it friendship.
// // My oh my, the devious crap we can pull with template traits!
// template<int nDim, typename TensorType> class GeomTensor;
// template<int nDim, typename TensorType> struct TensorFriendTraits {
//   typedef GeomTensor<nDim, TensorType> OtherTensorType; // Bogus!
// };
// template<>
// struct TensorFriendTraits<1, FullTensor> {
//   typedef GeomTensor<1, SymmetricTensor> OtherTensorType;
// };
// template<>
// struct TensorFriendTraits<2, FullTensor> {
//   typedef GeomTensor<2, SymmetricTensor> OtherTensorType;
// };
// template<>
// struct TensorFriendTraits<3, FullTensor> {
//   typedef GeomTensor<3, SymmetricTensor> OtherTensorType;
// };
// template<>
// struct TensorFriendTraits<1, SymmetricTensor> {
//   typedef GeomTensor<1, FullTensor> OtherTensorType;
// };
// template<>
// struct TensorFriendTraits<2, SymmetricTensor> {
//   typedef GeomTensor<2, FullTensor> OtherTensorType;
// };
// template<>
// struct TensorFriendTraits<3, SymmetricTensor> {
//   typedef GeomTensor<3, FullTensor> OtherTensorType;
// };

}

#endif

