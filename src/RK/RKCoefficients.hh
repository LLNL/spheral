//---------------------------------Spheral++----------------------------------//
// RKCoefficients
//
// Wraps the reproducing kernel correction coefficients
//----------------------------------------------------------------------------//
#ifndef __LLNLSpheral_RKCoefficients__
#define __LLNLSpheral_RKCoefficients__

#include "RK/RKCorrectionParams.hh"
#include <vector>
#include <iostream>

namespace Spheral {

template<typename Dimension>
struct RKCoefficients {
  RKOrder correctionOrder;        // The correction order
  std::vector<double> coeffs;     // The coefficients for use with ReproducingKernel/RKUtilities

  // Constructors and such
  RKCoefficients(): correctionOrder(RKOrder::ZerothOrder), coeffs() {}
  RKCoefficients(const RKCoefficients& rhs): correctionOrder(rhs.correctionOrder), coeffs(rhs.coeffs) {}
  RKCoefficients& operator=(const RKCoefficients& rhs) { correctionOrder = rhs.correctionOrder; coeffs = rhs.coeffs; return *this; }

  // Some convient methods to make us behave like a std::vector<double>
  typedef std::vector<double>::value_type value_type;
  typedef std::vector<double>::allocator_type allocator_type;
  typedef std::vector<double>::size_type size_type;
  typedef std::vector<double>::difference_type difference_type;
  typedef std::vector<double>::reference reference;
  typedef std::vector<double>::const_reference const_reference;
  typedef std::vector<double>::pointer pointer;
  typedef std::vector<double>::const_pointer const_pointer;
  typedef std::vector<double>::iterator iterator;
  typedef std::vector<double>::const_iterator const_iterator;
  typedef std::vector<double>::reverse_iterator reverse_iterator;
  typedef std::vector<double>::const_reverse_iterator const_reverse_iterator;
  
  reference at(size_type i)                                   { return coeffs.at(i); }    
  reference operator[](size_type i)                           { return coeffs[i]; }       
  reference front()                                           { return coeffs.front(); } 
  reference back()                                            { return coeffs.back(); }  
  pointer data()                                              { return coeffs.data(); }   
  iterator begin()                                            { return coeffs.begin(); }  
  iterator end()                                              { return coeffs.end(); }    
  reverse_iterator rbegin()                                   { return coeffs.rbegin(); } 
  reverse_iterator rend()                                     { return coeffs.rend(); }   

  const_reference at(size_type i)                       const { return coeffs.at(i); }
  const_reference operator[](size_type i)               const { return coeffs[i]; }
  const_reference front()                               const { return coeffs.front(); }
  const_reference back()                                const { return coeffs.back(); }
  const_pointer data()                                  const { return coeffs.data(); }
  const_iterator begin()                                const { return coeffs.begin(); }
  const_iterator end()                                  const { return coeffs.end(); }
  const_reverse_iterator rbegin()                       const { return coeffs.rbegin(); }
  const_reverse_iterator rend()                         const { return coeffs.rend(); }
  
  const_iterator cbegin()                               const { return coeffs.cbegin(); }
  const_iterator cend()                                 const { return coeffs.cend(); }
  const_reverse_iterator crbegin()                      const { return coeffs.crbegin(); }
  const_reverse_iterator crend()                        const { return coeffs.crend(); }

  void clear()                                                { coeffs.clear(); }
  bool empty()                                          const { return coeffs.empty(); }
  size_type size()                                      const { return coeffs.size(); }
  void resize(size_type count, double value=0.0)              { coeffs.resize(count, value); }
  
  bool operator==(const RKCoefficients<Dimension>& rhs) const { return correctionOrder == rhs.correctionOrder and coeffs == rhs.coeffs; }
  bool operator!=(const RKCoefficients<Dimension>& rhs) const { return correctionOrder != rhs.correctionOrder or  coeffs != rhs.coeffs; }
  bool operator< (const RKCoefficients<Dimension>& rhs) const { return (correctionOrder < rhs.correctionOrder ? true :
                                                                        correctionOrder > rhs.correctionOrder ? false :
                                                                        coeffs < rhs.coeffs); }
  bool operator> (const RKCoefficients<Dimension>& rhs) const { return (correctionOrder > rhs.correctionOrder ? true :
                                                                        correctionOrder < rhs.correctionOrder ? false :
                                                                        coeffs > rhs.coeffs); }
  bool operator<=(const RKCoefficients<Dimension>& rhs) const { return (*this) == rhs or (*this) < rhs; }
  bool operator>=(const RKCoefficients<Dimension>& rhs) const { return (*this) == rhs or (*this) > rhs; }
};

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::ostream&
operator<<(std::ostream& os, const RKCoefficients<Dimension>& x) {
  os << "[";
  if (not x.empty()) {
    const auto n1 = x.size() - 1;
    for (auto i = 0; i < n1; ++i) os << x[i] << " ";
    os << x.back();
  }
  os << "]";
  return os;
}

}

#else

namespace Spheral {
  template<typename Dimension> struct RKCoefficients;
}

#endif

  
