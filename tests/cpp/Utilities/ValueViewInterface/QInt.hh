//-----------------------------------------------------------------------------
//  Quadratic Interpolator Example implementation
//-----------------------------------------------------------------------------

#include "Utilities/ValueViewInterface.hh"

namespace Spheral {

// The Data class needs to be CHAICopyable in order to trigger nested copies for Copyable 
// members within.
class QIntData : public SPHERALCopyable<QIntData>{
public:

  SPHERAL_HOST_DEVICE QIntData() = default;
  SPHERAL_HOST_DEVICE QIntData(QIntData const& rhs) = default;
  SPHERAL_HOST_DEVICE QIntData& operator=(QIntData const& rhs) = default;

  using CoeffsType = ManagedVector<double>;

  double mXmin, mXmax, mXstep;
  CoeffsType mcoeffs;

  SPHERAL_HOST void initialize(size_t min)
  {
    mXmin = min;
    mcoeffs.resize(10);
    mcoeffs[0] = 0.1;
    mcoeffs[1] = 0.2;
    mcoeffs[2] = 0.3;
    mcoeffs[9] = 0.19;
  }

  SPHERAL_HOST void editData(size_t min)
  {
    mXmin = min;
    mcoeffs[9] = 91;
  }

  SPHERAL_HOST_DEVICE double xmin() const { return mXmin; }
  SPHERAL_HOST_DEVICE double xmax() const { return mXmax; }
  SPHERAL_HOST_DEVICE CoeffsType const& coeffs() const { return mcoeffs; }

  // Define the required interface for a SPHERALCopyable object.
  void free() { mcoeffs.free(); }
  SPHERAL_HOST_DEVICE QIntData& operator=(std::nullptr_t) { mcoeffs=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(QIntData const& rhs) { *this = rhs; }
};


class QInt;

class QIntView : public SpheralViewInterface<QIntData>
{
  VIEW_INTERFACE_METACLASS((QInt), (QIntView), (QIntData))
public:
  friend class QInt;
  using CoeffsType = typename QIntData::CoeffsType;
protected:
  //VIEW_DEFINE_ALLOC_CTOR(QIntView)
  // Interal interface for accessing the underlying members of QIntData
  SMART_PTR_MEMBER_ACCESSOR(CoeffsType, mcoeffs)

public:
  // Forward View capable methods
  SPHERAL_HOST_DEVICE double xmin() const { return sptr_data().xmin(); }
  SPHERAL_HOST_DEVICE double xmax() const { return sptr_data().xmax(); }
  SPHERAL_HOST_DEVICE CoeffsType const& coeffs() const { return sptr_data().coeffs(); }
};



class QInt : public SpheralValueInterface<QIntView>
{
  VALUE_TYPE_ALIASES((QIntView))
public:
  VALUE_DEF_CTOR(QInt)
  //VALUE_COPY_CTOR(QInt, QIntData)
  //VALUE_ASSIGNEMT_OP(QInt, QIntData)
  VALUE_TOVIEW_OP()

  // Forward Value capable methods
  SPHERAL_HOST void initialize(size_t min) const { return sptr_data().initialize(min); }
  SPHERAL_HOST void editData(size_t min) const { return sptr_data().editData(min); }
  SPHERAL_HOST CoeffsType coeffs() const { return deepCopy(sptr_data().coeffs()); }
};

} // namespace Spheral

