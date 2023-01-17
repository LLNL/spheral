#ifndef VEC_HH
#define VEC_HH

#if USE_DEVICE
#define CAST_ATOMIC_TYPE(X) reinterpret_cast<camp::decay<decltype(X)>::ATOMIC_TYPE&>(X) 
#else
#define CAST_ATOMIC_TYPE(X) X
#endif

class vec_base {
public:
  using DATA_TYPE = double;

  RAJA_HOST_DEVICE vec_base(const DATA_TYPE d) :
    mx(d)
#if VEC3
    ,my(d)
    ,mz(d)
#endif
  {};

protected:
    DATA_TYPE mx = 0;
#if VEC3
    DATA_TYPE my = 0;
    DATA_TYPE mz = 0;
#endif
};

class atomic_vec;

class vec : public vec_base{
public:

  using ATOMIC_TYPE = atomic_vec;

  RAJA_HOST_DEVICE vec(const DATA_TYPE d = 0) : vec_base(d) {};
  RAJA_HOST_DEVICE vec(const vec& pd) : vec_base(pd.mx) {};
  

  std::string string() const {return std::to_string(mx);}

  RAJA_HOST_DEVICE DATA_TYPE get_data() const {return mx;}


  RAJA_HOST_DEVICE inline vec& operator+=(const vec& rhs) {
    this->mx += rhs.mx;
#if VEC3
    this->my += rhs.my;
    this->mz += rhs.mz;
#endif
    return *this;
  }
  RAJA_HOST_DEVICE inline vec operator+(const vec& rhs) const {
    vec result(*this);
    result+=rhs;
    return result;
  }

  RAJA_HOST_DEVICE bool operator!=(const vec& rhs) const {return !(mx == rhs.mx);}

};

class atomic_vec : public vec{
public:
  RAJA_HOST_DEVICE atomic_vec(const DATA_TYPE d = 0) : vec(d) {};
  RAJA_HOST_DEVICE atomic_vec(const vec& pd) : vec(pd) {};

  RAJA_HOST_DEVICE inline atomic_vec operator+=(const atomic_vec& rhs) {
    ATOMIC_ADD(&mx, rhs.mx);
#if VEC3
    ATOMIC_ADD(&my, rhs.my);
    ATOMIC_ADD(&mz, rhs.mz);
#endif
    return *this;
  }

  RAJA_HOST_DEVICE inline atomic_vec operator+(const atomic_vec& rhs) {
    atomic_vec result(*this);
    result+=rhs;
    return result;
  }

};

std::ostream& operator<<(std::ostream& out, const vec& d) {
  return out << d.string();
}

#endif
