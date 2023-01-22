//---------------------------------Spheral++----------------------------------//
// ManufacturedSolution
//
// Represents some manufactured solutions along with integration coefficients
//----------------------------------------------------------------------------//
#ifndef __Spheral_ManufacturedSolution_hh__
#define __Spheral_ManufacturedSolution_hh__

#include "IntegrationCoefficient.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// ManufacturedFunction
//
// Lets us specify a given function along with its derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
class ManufacturedFunction : public IntegrationCoefficient<Dimension, double> {
public:
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  
  ManufacturedFunction();
  virtual double evaluate(const double t, const Vector& x) const = 0;
  virtual Vector evaluateSpatialGradient(const double t, const Vector& x) const = 0;
  virtual SymTensor evaluateSpatialHessian(const double t, const Vector& x) const = 0;
  virtual double evaluateTimeDerivative(const double t, const Vector& x) const = 0;
  virtual double evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const override;
};

//------------------------------------------------------------------------------
// ManufacturedSteadyStateFunction
//
// k_{0}
//------------------------------------------------------------------------------
template<typename Dimension>
class ManufacturedSteadyStateFunction : public ManufacturedFunction<Dimension> {
public:
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  
  ManufacturedSteadyStateFunction(const double t,
                                  std::shared_ptr<ManufacturedFunction<Dimension>> func);
  virtual double evaluate(const double t, const Vector& x) const override;
  virtual Vector evaluateSpatialGradient(const double t, const Vector& x) const override;
  virtual SymTensor evaluateSpatialHessian(const double t, const Vector& x) const override;
  virtual double evaluateTimeDerivative(const double t, const Vector& x) const override;
private:
  double mT;
  std::shared_ptr<ManufacturedFunction<Dimension>> mFunc;
};

//------------------------------------------------------------------------------
// ManufacturedConstantFunction
//
// k_{0}
//------------------------------------------------------------------------------
template<typename Dimension>
class ManufacturedConstantFunction : public ManufacturedFunction<Dimension> {
public:
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  
  ManufacturedConstantFunction(const double coefficient);
  virtual double evaluate(const double t, const Vector& x) const override;
  virtual Vector evaluateSpatialGradient(const double t, const Vector& x) const override;
  virtual SymTensor evaluateSpatialHessian(const double t, const Vector& x) const override;
  virtual double evaluateTimeDerivative(const double t, const Vector& x) const override;
private:
  const double mCoefficient;
};

//------------------------------------------------------------------------------
// ManufacturedSinusoidalFunction
//
// k_{0}\left(1+\frac{1}{2\pi}\prod_{\alpha}\cos\left(k_{\alpha}\left(x^{\alpha}+t\right)\right)\right)
//------------------------------------------------------------------------------
template<typename Dimension>
class ManufacturedSinusoidalFunction : public ManufacturedFunction<Dimension> {
public:
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  
  ManufacturedSinusoidalFunction(const std::vector<double>& coefficients);
  virtual double evaluate(const double t, const Vector& x) const override;
  virtual Vector evaluateSpatialGradient(const double t, const Vector& x) const override;
  virtual SymTensor evaluateSpatialHessian(const double t, const Vector& x) const override;
  virtual double evaluateTimeDerivative(const double t, const Vector& x) const override;
private:
  std::vector<double> mCoefficients;
};
  
//------------------------------------------------------------------------------
// ManufacturedWaveFunction
//
// k_{0}\left[1+\frac{1}{t^{2}+6}\exp\left(-k_{1}\left[\left|x\right|-t\right]^{2}\right)\right]
//------------------------------------------------------------------------------
template<typename Dimension>
class ManufacturedWaveFunction : public ManufacturedFunction<Dimension> {
public:
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  
  ManufacturedWaveFunction(const std::vector<double>& coefficients);
  virtual double evaluate(const double t, const Vector& x) const override;
  virtual Vector evaluateSpatialGradient(const double t, const Vector& x) const override;
  virtual SymTensor evaluateSpatialHessian(const double t, const Vector& x) const override;
  virtual double evaluateTimeDerivative(const double t, const Vector& x) const override;
private:
  std::vector<double> mCoefficients;
};
  
//------------------------------------------------------------------------------
// ManufacturedTransportSolution
//
// Calculates a source for the given coefficients 
//------------------------------------------------------------------------------
template<typename Dimension>
class ManufacturedTransportSolution : public IntegrationCoefficient<Dimension, std::vector<double>> {
public:
  typedef typename Dimension::Vector Vector;
  
  ManufacturedTransportSolution(const double c,
                                const int numOrdinates,
                                const double angularNorm,
                                const std::vector<Vector>& ordinates,
                                std::shared_ptr<ManufacturedFunction<Dimension>> psiFunc,
                                std::shared_ptr<ManufacturedFunction<Dimension>> sigmaAFunc);
  virtual double evaluatePhi(const double t, const Vector& x) const;
  virtual std::vector<double> evaluatePsi(const double t, const Vector& x) const;
  virtual std::vector<double> evaluateSource(const double t, const Vector& x) const;
  virtual std::vector<double> evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const override;
  
private:
  double mCInv;
  int mNumOrdinates;
  double mAngularNorm;
  const std::vector<Vector>& mOrdinates;
  std::shared_ptr<ManufacturedFunction<Dimension>> mPsiFunc;
  std::shared_ptr<ManufacturedFunction<Dimension>> mSigmaAFunc;
  std::shared_ptr<ManufacturedFunction<Dimension>> mSigmaSFunc;
};

} // end namespace Spheral

#endif
