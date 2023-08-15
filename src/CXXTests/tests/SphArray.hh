#ifndef SPHARRAY_H
#define SPHARRAY_H

#include "LvArray/ChaiBuffer.hpp"

class ArrayContainerType {
  using DataType = double;
  using ContainerType = chai::ManagedArray<DataType>;
  ContainerType mdata;
public:
  ArrayContainerType() : mdata() {}

  ArrayContainerType(const std::vector<double>& vals) {
    mdata.allocate(vals.size(), chai::CPU);
    mdata.registerTouch(chai::CPU);
    for (int i = 0 ; i < vals.size(); i++) mdata[i] = vals[i];
  }

  RAJA_HOST_DEVICE
  ArrayContainerType(const ArrayContainerType& rhs) : mdata(rhs.mdata) {}

  RAJA_HOST_DEVICE
  const DataType& operator[](size_t i) const { return mdata[i]; }

  RAJA_HOST_DEVICE
  size_t size() const { return mdata.size(); }

  RAJA_HOST_DEVICE
  const ContainerType& coeffs() const {return mdata;}
  
private:
};


#endif // SPHARRAY_H
