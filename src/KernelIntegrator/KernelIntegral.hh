//---------------------------------Spheral++----------------------------------//
// KernelIntegral
//
// Represents an integral to be integrated by KernelIntegrator
//----------------------------------------------------------------------------//
#ifndef __Spheral_KernelIntegral_hh__
#define __Spheral_KernelIntegral_hh__

#define REPLACEOVERLAP true

#include <memory>
#include <vector>
#include <unordered_map>
#include "FlatConnectivity.hh"
#include "IntegrationCoefficient.hh"
#include "KernelIntegrationData.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// KernelIntegralBase
// 
// Base integral class
//------------------------------------------------------------------------------
template<typename Dimension>
class KernelIntegralBase {
public:
  KernelIntegralBase() { }

  // Does integral depend on bilinear indexing?
  virtual bool bilinear() const = 0;

  // Does integral have volume integral terms?
  virtual bool volume() const = 0;

  // Does integral have surface integral terms?
  virtual bool surface() const = 0;
  
  // Zero out the data members and resize them if needed
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) = 0;

  // Perform any post-integration operations
  virtual void finalize(const FlatConnectivity<Dimension>& flatConnectivity) { }
  
  // Add a value for the given integration point to the volume integral
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) { }

  // Add a value for the given integration point to the surface integral
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) { }
};

//------------------------------------------------------------------------------
// KernelIntegral
// 
// Includes return type for data
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
class KernelIntegral : public KernelIntegralBase<Dimension> {
public:
  typedef std::vector<DataType> StorageType;
  
  KernelIntegral() { }

  // Return the value of the integral
  virtual const StorageType& values() const = 0;
};

//------------------------------------------------------------------------------
// LinearIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
class LinearIntegral : public KernelIntegral<Dimension, DataType> {
public:
  typedef std::vector<DataType> StorageType;
  
  LinearIntegral() { }

  // Inherited from KernelIntegral
  virtual bool bilinear() const override { return false; }
  virtual const StorageType& values() const override { return mValues; }
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) override;
  
protected:
  bool mVolume;
  bool mSurface;
  StorageType mValues;
};

//------------------------------------------------------------------------------
// BilinearIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
class BilinearIntegral : public KernelIntegral<Dimension, std::vector<BaseDataType>> {
public:
  typedef std::vector<std::vector<BaseDataType>> StorageType;
  
  BilinearIntegral() { }
  
  // Inherited from KernelIntegral
  virtual bool bilinear() const override { return true; }
  virtual const StorageType& values() const override { return mValues; }
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) override;
  
protected:
  StorageType mValues;
};

//------------------------------------------------------------------------------
// LinearSurfaceDependentIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
class LinearSurfaceDependentIntegral : public KernelIntegral<Dimension, std::vector<BaseDataType>> {
public:
  typedef std::vector<std::vector<BaseDataType>> StorageType;
  
  LinearSurfaceDependentIntegral() { }

  // Inherited from KernelIntegral
  virtual bool bilinear() const override { return false; }
  virtual bool volume() const override { return false; }
  virtual bool surface() const override { return true; }
  virtual const StorageType& values() const override { return mValues; }
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) override;
  
protected:
  StorageType mValues;
};

//------------------------------------------------------------------------------
// BilinearSurfaceDependentIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
class BilinearSurfaceDependentIntegral : public KernelIntegral<Dimension, std::vector<BaseDataType>> {
public:
  typedef std::vector<std::vector<BaseDataType>> StorageType;
  
  BilinearSurfaceDependentIntegral() { }
  
  virtual bool bilinear() const override { return true; }
  virtual bool volume() const override { return false; }
  virtual bool surface() const override { return true; }
  virtual const StorageType& values() const override { return mValues; }
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) override;
  
protected:
  StorageType mValues;
};

//------------------------------------------------------------------------------
// LinearKernel
// 
// Integral of kernel with a functional coefficient
// \int_{V}ku_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearKernel : public LinearIntegral<Dimension, typename Dimension::Scalar>,
                     public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  LinearKernel() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// LinearGrad
// 
// Integral of kernel with a functional coefficient
// \int_{V}k\partial^{\alpha}u_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearGrad : public LinearIntegral<Dimension, typename Dimension::Vector>,
                   public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  LinearGrad() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// LinearKernelVector
// 
// Integral of kernel with a vector coefficient
// \int_{V}k^{\alpha}u_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearKernelVector : public LinearIntegral<Dimension, typename Dimension::Vector>,
                           public IntegralDependsOnCoefficient<Dimension, typename Dimension::Vector> {
public:
  LinearKernelVector() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// LinearKernelStdVector
// 
// Integral of kernel with a functional coefficient
// \int_{V}ku_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearKernelStdVector : public LinearIntegral<Dimension, std::vector<typename Dimension::Scalar>>,
                              public IntegralDependsOnCoefficient<Dimension, std::vector<typename Dimension::Scalar>> {
public:
  LinearKernelStdVector(size_t size): mSize(size) { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) override;

private:
  size_t mSize;
};

//------------------------------------------------------------------------------
// LinearGradStdVector
// 
// Integral of kernel with a functional coefficient
// \int_{V}k\partial^{\alpha}u_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearGradStdVector : public LinearIntegral<Dimension, std::vector<typename Dimension::Vector>>,
                            public IntegralDependsOnCoefficient<Dimension, std::vector<typename Dimension::Scalar>> {
public:
  LinearGradStdVector(size_t size): mSize(size) { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) override;

private:
  size_t mSize;
};

//------------------------------------------------------------------------------
// BilinearKernelKernel
// 
// Bilinear integral of kernel and kernel
// \int_{V}ku_{i}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearKernelKernel : public BilinearIntegral<Dimension, typename Dimension::Scalar>,
                             public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearKernelKernel() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearGradKernel
// 
// Bilinear integral of grad kernel and kernel
// \int_{V}k\partial^{\alpha}u_{i}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearGradKernel : public BilinearIntegral<Dimension, typename Dimension::Vector>,
                           public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearGradKernel() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearKernelGrad
// 
// Bilinear integral of grad kernel and kernel
// \int_{V}ku_{i}\partial^{\alpha}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearKernelGrad : public BilinearIntegral<Dimension, typename Dimension::Vector>,
                           public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearKernelGrad() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearGradDotGrad
// 
// Bilinear integral of grad kernel and grad kernel
// \int_{V}k\partial^{\alpha}u_{i}\partial^{\alpha}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearGradDotGrad : public BilinearIntegral<Dimension, typename Dimension::Scalar>,
                            public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearGradDotGrad() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearGradProdGrad
// 
// Bilinear integral of grad kernel and grad kernel
// \int_{V}k\partial^{\alpha}u_{i}\partial^{\beta}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearGradProdGrad : public BilinearIntegral<Dimension, typename Dimension::Tensor>,
                             public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearGradProdGrad() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearSurfaceNormalKernelKernelFromGrad
// 
// Bilinear surface integral of kernel and kernel times normal, calculated using gradients
// \int_{S}n^{\alpha}u_{i}u_{j}=\int_{V}\partial^{\alpha}u_{i}u_{j}+\int_{V}u_{i}\partial^{\alpha}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearSurfaceNormalKernelKernelFromGrad : public BilinearIntegral<Dimension, typename Dimension::Vector>,
                                                  public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearSurfaceNormalKernelKernelFromGrad() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// LinearSurfaceKernel
// 
// Integral of kernel over a surface with a functional coefficient
// \int_{S}ku_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearSurfaceKernel : public LinearIntegral<Dimension, typename Dimension::Scalar>,
                            public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  LinearSurfaceKernel() { }
  
  virtual bool volume() const override { return false; }
  virtual bool surface() const override { return true; }
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// LinearSurfaceNormalKernel
// 
// Linear surface integral of kernel and normal
// \int_{S}kn^{\alpha}u_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearSurfaceNormalKernel : public LinearSurfaceDependentIntegral<Dimension, typename Dimension::Vector>,
                                  public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  LinearSurfaceNormalKernel() { }
  
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// LinearSurfaceNormalKernelStdVector
// 
// Linear surface integral of kernel and normal
// \int_{S}kn^{\alpha}u_{i}
//------------------------------------------------------------------------------
template<typename Dimension>
class LinearSurfaceNormalKernelStdVector : public LinearSurfaceDependentIntegral<Dimension, std::vector<typename Dimension::Vector>>,
                                           public IntegralDependsOnCoefficient<Dimension, std::vector<typename Dimension::Scalar>> {
public:
  LinearSurfaceNormalKernelStdVector(size_t size): mSize(size) { }
  
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) override;
  virtual void initialize(const FlatConnectivity<Dimension>& flatConnectivity) override;

private:
  size_t mSize;
};

//------------------------------------------------------------------------------
// BilinearSurfaceKernelKernel
// 
// Bilinear surface integral of kernel and kernel
// \int_{S}ku_{i}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearSurfaceKernelKernel : public BilinearIntegral<Dimension, typename Dimension::Scalar>,
                                    public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearSurfaceKernelKernel() { }
  
  virtual bool volume() const override { return false; }
  virtual bool surface() const override { return true; }
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearSurfaceNormalKernelKernel
// 
// Bilinear surface integral of kernel, kernel, and normal
// \int_{S}kn^{\alpha}u_{i}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearSurfaceNormalKernelKernel : public BilinearSurfaceDependentIntegral<Dimension, typename Dimension::Vector>,
                                          public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearSurfaceNormalKernelKernel() { }
  
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearSurfaceNormalKernelDotGrad
// 
// Bilinear surface integral of kernel and normal dotted into grad
// \int_{S}kn^{\alpha}u_{i}\partial^{\alpha}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension>
class BilinearSurfaceNormalKernelDotGrad : public BilinearIntegral<Dimension, typename Dimension::Scalar>,
                                           public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearSurfaceNormalKernelDotGrad() { }
  
  virtual bool volume() const override { return false; }
  virtual bool surface() const override { return true; }
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// CellCoefficient
// 
// Integral of an arbitrary coefficient in each cell
// \int_{V}k
//------------------------------------------------------------------------------
template<typename Dimension>
class CellCoefficient : public LinearIntegral<Dimension, typename Dimension::Scalar>,
                        public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  CellCoefficient() { }
  
  virtual bool volume() const override { return true; }
  virtual bool surface() const override { return false; }
  virtual void addToIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// SurfaceNormalCoefficient
// 
// Integral of an arbitrary coefficient over each surface
// \int_{S}n^{\alpha}k
// Since our storage is sized for the overlap of a kernel with surfaces, we could
// make this more efficient if it were ever used outside of testing
//------------------------------------------------------------------------------
template<typename Dimension>
class SurfaceNormalCoefficient : public LinearSurfaceDependentIntegral<Dimension, typename Dimension::Vector>,
                                 public IntegralDependsOnCoefficient<Dimension, typename Dimension::Scalar> {
public:
  SurfaceNormalCoefficient() { }
  
  virtual void addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) override;
};

//------------------------------------------------------------------------------
// BilinearMultiplyByFieldList
// 
// Bilinear integral, multiplied by a FieldList
// If the original integral is
// \int_{V}u_{i}u_{j},
// then this returns
// f_{i}\int_{V}u_{i}u_{j}
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
class BilinearMultiplyByFieldList : public BilinearIntegral<Dimension, BaseDataType>,
                                    public IntegralDependsOnFieldListCoefficient<Dimension, typename Dimension::Scalar> {
public:
  BilinearMultiplyByFieldList(std::shared_ptr<BilinearIntegral<Dimension, BaseDataType>> baseIntegral) : mBaseIntegral(baseIntegral) { }
  
  virtual bool volume() const override { return mBaseIntegral->volume(); }
  virtual bool surface() const override { return mBaseIntegral->surface(); }
  virtual void finalize(const FlatConnectivity<Dimension>& flatConnectivity) override;
private:
  std::shared_ptr<BilinearIntegral<Dimension, BaseDataType>> mBaseIntegral;
};

} // end namespace Spheral

#include "KernelIntegralInline.hh"

#endif
