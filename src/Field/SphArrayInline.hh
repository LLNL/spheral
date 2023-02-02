namespace Spheral {

template<typename sph_array_t>
inline
SphArrayIterator<sph_array_t>::
SphArrayIterator(pointer ptr)
  : mPtr(ptr) {}

template<typename sph_array_t>
inline
SphArrayIterator<sph_array_t>&
SphArrayIterator<sph_array_t>::
operator++() { mPtr++; return *this; }

template<typename sph_array_t>
inline
SphArrayIterator<sph_array_t>
SphArrayIterator<sph_array_t>::
operator++(int) { SphArrayIterator it = *this; operator++(); return it; }

template<typename sph_array_t>
inline
SphArrayIterator<sph_array_t>&
SphArrayIterator<sph_array_t>::
operator--() { mPtr--; return *this; }

template<typename sph_array_t>
inline
SphArrayIterator<sph_array_t>
SphArrayIterator<sph_array_t>::
operator--(int) { SphArrayIterator it = *this; operator--(); return it; }

template<typename sph_array_t>
inline
SphArrayIterator<sph_array_t>
SphArrayIterator<sph_array_t>::
operator+(const difference_type& index) const { SphArrayIterator it = *this; it.mPtr += index; return it; }

template<typename sph_array_t>
inline
SphArrayIterator<sph_array_t>
SphArrayIterator<sph_array_t>::
operator-(const difference_type& index) const { SphArrayIterator it = *this; it.mPtr -= index; return it; }

template<typename sph_array_t>
inline
typename SphArrayIterator<sph_array_t>::difference_type
SphArrayIterator<sph_array_t>::
operator-(const SphArrayIterator& it){ return std::distance(it.mPtr,mPtr); }

template<typename sph_array_t>
inline
typename SphArrayIterator<sph_array_t>::reference
SphArrayIterator<sph_array_t>::
operator[](int index){ return *(mPtr+index); }

template<typename sph_array_t>
inline
typename SphArrayIterator<sph_array_t>::pointer
SphArrayIterator<sph_array_t>::
operator->(){ return mPtr; }

template<typename sph_array_t>
inline
typename SphArrayIterator<sph_array_t>::reference
SphArrayIterator<sph_array_t>::
operator*() const { return *mPtr; }

template<typename sph_array_t>
inline
bool
SphArrayIterator<sph_array_t>::
operator==(const SphArrayIterator& rhs) const { return mPtr == rhs.mPtr; }

template<typename sph_array_t>
inline
bool
SphArrayIterator<sph_array_t>::
operator!=(const SphArrayIterator& rhs) const { return !(*this == rhs); }

template<typename sph_array_t>
inline
bool
SphArrayIterator<sph_array_t>::
operator<(const SphArrayIterator& rhs) const {return mPtr < rhs.mPtr; }

template<typename sph_array_t>
inline
bool
SphArrayIterator<sph_array_t>::
operator<=(const SphArrayIterator& rhs) const {return operator<(rhs) || operator==(rhs); }

template<typename sph_array_t>
inline
bool
SphArrayIterator<sph_array_t>::
operator>(const SphArrayIterator& rhs) const {return !(operator<(rhs) || operator==(rhs)); }

template<typename sph_array_t>
inline
bool
SphArrayIterator<sph_array_t>::
operator>=(const SphArrayIterator& rhs) const {return !operator<(rhs); }

template<typename sph_array_t>
SphArrayIterator<sph_array_t>::
operator SphArrayIterator<typename sph_array_t::ViewTypeConst>() const {return SphArrayIterator<typename sph_array_t::ViewTypeConst>(mPtr); }

} //  namespace Spheral
