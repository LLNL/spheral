#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "Field/Field.hh"
#include "Geometry/Dimension.hh"
#include "DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Generic Field -> string conversion.
//------------------------------------------------------------------------------
template<typename D, typename T>
inline
std::string
field2string(const FieldSpace::Field<D, T>& val,
             const int precision) {
  return val.string(precision);
}

//------------------------------------------------------------------------------
// Generic string -> Field conversion.
//------------------------------------------------------------------------------
template<typename D, typename T>
inline
void
string2field(FieldSpace::Field<D, T>& field,
             const std::string& val) {
  field.string(val);
}

//------------------------------------------------------------------------------
// Specializations: Field<D, int> -> string.
//------------------------------------------------------------------------------
template<>
inline
std::string
field2string<Dim<1>, int>(const FieldSpace::Field<Dim<1>, int>& val,
                          const int precision) {
  std::string result;

  char tmp[21];
  tmp[20] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%20i ", val[i]);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<2>, int>(const FieldSpace::Field<Dim<2>, int>& val,
                          const int precision) {
  std::string result;

  char tmp[21];
  tmp[20] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%20i ", val[i]);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<3>, int>(const FieldSpace::Field<Dim<3>, int>& val,
                          const int precision) {
  std::string result;

  char tmp[21];
  tmp[20] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%20i ", val[i]);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// Specializations: Field<D, Scalar> -> string.
//------------------------------------------------------------------------------
template<>
inline
std::string
field2string<Dim<1>, Dim<1>::Scalar>(const FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& val,
                                     const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e ", val[i]);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<2>, Dim<2>::Scalar>(const FieldSpace::Field<Dim<2>, Dim<2>::Scalar>& val,
                                     const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e ", val[i]);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<3>, Dim<3>::Scalar>(const FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& val,
                                     const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e ", val[i]);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// Specializations: Field<D, Vector> -> string.
//------------------------------------------------------------------------------
template<>
inline
std::string
field2string<Dim<1>, Dim<1>::Vector>(const FieldSpace::Field<Dim<1>, Dim<1>::Vector>& val,
                                     const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e ", val[i].x());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<2>, Dim<2>::Vector>(const FieldSpace::Field<Dim<2>, Dim<2>::Vector>& val,
                                     const int precision) {
  std::string result;

  char tmp[31*2 + 1];
  tmp[31*2] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e %30.20e ", val[i].x(), val[i].y());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<3>, Dim<3>::Vector>(const FieldSpace::Field<Dim<3>, Dim<3>::Vector>& val,
                                     const int precision) {
  std::string result;

  char tmp[31*3 + 1];
  tmp[31*3] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e %30.20e %30.20e ", val[i].x(), val[i].y(), val[i].z());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// Specializations: Field<D, Tensor> -> string.
//------------------------------------------------------------------------------
template<>
inline
std::string
field2string<Dim<1>, Dim<1>::Tensor>(const FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& val,
                                     const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e ", val[i].xx());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<2>, Dim<2>::Tensor>(const FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& val,
                                     const int precision) {
  std::string result;

  char tmp[31*4 + 1];
  tmp[31*4] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e %30.20e %30.20e %30.20e ", 
            val[i].xx(), val[i].xy(),
            val[i].yx(), val[i].yy());
    result += std::string(tmp);
  }

  return result;
}

template<>
inline
std::string
field2string<Dim<3>, Dim<3>::Tensor>(const FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& val,
                                     const int precision) {
  std::string result;

  char tmp[31*9 + 1];
  tmp[31*9] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e ", 
            val[i].xx(), val[i].xy(), val[i].xz(), 
            val[i].yx(), val[i].yy(), val[i].yz(), 
            val[i].zx(), val[i].zy(), val[i].zz());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// Specializations: Field<D, SymTensor> -> string.
//------------------------------------------------------------------------------
template<>
inline
std::string
field2string<Dim<1>, Dim<1>::SymTensor>(const FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& val,
                                        const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e ", val[i].xx());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<2>, Dim<2>::SymTensor>(const FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& val,
                                        const int precision) {
  std::string result;

  char tmp[31*3 + 1];
  tmp[31*3] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e %30.20e %30.20e ", 
            val[i].xx(), val[i].xy(),
                         val[i].yy());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

template<>
inline
std::string
field2string<Dim<3>, Dim<3>::SymTensor>(const FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& val,
                                        const int precision) {
  std::string result;

  char tmp[31*6 + 1];
  tmp[31*6] = 0;
  int n = val.nodeList().numInternalNodes();
  sprintf(tmp, "%20i", n);
  result += std::string(tmp);

  for (int i = 0; i != n; ++i) {
    sprintf(tmp, "%30.20e %30.20e %30.20e %30.20e %30.20e %30.20e ", 
            val[i].xx(), val[i].xy(), val[i].xz(), 
                         val[i].yy(), val[i].yz(), 
                                      val[i].zz());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// Specializations: Field<D, ThirdRankTensor> -> string.
//------------------------------------------------------------------------------
template<>
inline
std::string
field2string<Dim<1>, Dim<1>::ThirdRankTensor>(const FieldSpace::Field<Dim<1>, Dim<1>::ThirdRankTensor>& val,
                                              const int precision) {
  typedef Spheral::Dim<1>::ThirdRankTensor Element;
  std::stringstream result;

  // Output the size of the field.
  const int n = val.nodeList().numInternalNodes();
  result << n << " ";

  // Now all the indvidual elements.
  result << std::setprecision(20) << std::scientific;
  for (int i = 0; i != n; ++i) result << val(i) << " ";

  return result.str();
}

template<>
inline
std::string
field2string<Dim<2>, Dim<2>::ThirdRankTensor>(const FieldSpace::Field<Dim<2>, Dim<2>::ThirdRankTensor>& val,
                                              const int precision) {
  typedef Spheral::Dim<2>::ThirdRankTensor Element;
  std::stringstream result;

  // Output the size of the field.
  const int n = val.nodeList().numInternalNodes();
  result << n << " ";

  // Now all the indvidual elements.
  result << std::setprecision(20) << std::scientific;
  for (int i = 0; i != n; ++i) result << val(i) << " ";

  return result.str();
}

template<>
inline
std::string
field2string<Dim<3>, Dim<3>::ThirdRankTensor>(const FieldSpace::Field<Dim<3>, Dim<3>::ThirdRankTensor>& val,
                                              const int precision) {
  typedef Spheral::Dim<3>::ThirdRankTensor Element;
  std::stringstream result;

  // Output the size of the field.
  const int n = val.nodeList().numInternalNodes();
  result << n << " ";

  // Now all the indvidual elements.
  result << std::setprecision(20) << std::scientific;
  for (int i = 0; i != n; ++i) result << val(i) << " ";

  return result.str();
}

//------------------------------------------------------------------------------
// Specializations: string -> Field<D, int>.
//------------------------------------------------------------------------------
template<>
inline
void
string2field<Dim<1>, int>(FieldSpace::Field<Dim<1>, int>& field,
                          const std::string& val) {
  const int len = 21;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    sscanf(valptr, "%i", &field(i));
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<2>, int>(FieldSpace::Field<Dim<2>, int>& field,
                          const std::string& val) {
  const int len = 21;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    sscanf(valptr, "%i", &field(i));
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<3>, int>(FieldSpace::Field<Dim<3>, int>& field,
                          const std::string& val) {
  const int len = 21;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    sscanf(valptr, "%i", &field(i));
  }
  delete [] valptr;
}

//------------------------------------------------------------------------------
// Specializations: string -> Field<D, Scalar>.
//------------------------------------------------------------------------------
template<>
inline
void
string2field<Dim<1>, Dim<1>::Scalar>(FieldSpace::Field<Dim<1>, Dim<1>::Scalar>& field,
                                     const std::string& val) {
  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    sscanf(valptr, "%lf", &field(i));
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<2>, Dim<2>::Scalar>(FieldSpace::Field<Dim<2>, Dim<2>::Scalar>& field,
                                     const std::string& val) {
  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    sscanf(valptr, "%lf", &field(i));
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<3>, Dim<3>::Scalar>(FieldSpace::Field<Dim<3>, Dim<3>::Scalar>& field,
                                     const std::string& val) {
  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    sscanf(valptr, "%lf", &field(i));
  }
  delete [] valptr;
}

//------------------------------------------------------------------------------
// Specializations: string -> Field<D, Vector>.
//------------------------------------------------------------------------------
template<>
inline
void
string2field<Dim<1>, Dim<1>::Vector>(FieldSpace::Field<Dim<1>, Dim<1>::Vector>& field,
                                     const std::string& val) {
  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double x;
    sscanf(valptr, "%lf", &x);
    field(i).x(x);
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<2>, Dim<2>::Vector>(FieldSpace::Field<Dim<2>, Dim<2>::Vector>& field,
                                     const std::string& val) {
  const int len = 31*2;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double x, y;
    sscanf(valptr, "%lf %lf", &x, &y);
    field(i).x(x);
    field(i).y(y);
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<3>, Dim<3>::Vector>(FieldSpace::Field<Dim<3>, Dim<3>::Vector>& field,
                                     const std::string& val) {
  const int len = 31*3;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double x, y, z;
    sscanf(valptr, "%lf %lf %lf", &x, &y, &z);
    field(i).x(x);
    field(i).y(y);
    field(i).z(z);
  }
  delete [] valptr;
}

//------------------------------------------------------------------------------
// Specializations: string -> Field<D, Tensor>.
//------------------------------------------------------------------------------
template<>
inline
void
string2field<Dim<1>, Dim<1>::Tensor>(FieldSpace::Field<Dim<1>, Dim<1>::Tensor>& field,
                                     const std::string& val) {
  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx;
    sscanf(valptr, "%lf", &xx);
    field(i).xx(xx);
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<2>, Dim<2>::Tensor>(FieldSpace::Field<Dim<2>, Dim<2>::Tensor>& field,
                                     const std::string& val) {
  const int len = 31*4;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, yx, yy;
    sscanf(valptr, "%lf %lf %lf %lf", &xx, &xy, &yx, &yy);
    field(i).xx(xx);
    field(i).xy(xy);
    field(i).yx(yx);
    field(i).yy(yy);
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<3>, Dim<3>::Tensor>(FieldSpace::Field<Dim<3>, Dim<3>::Tensor>& field,
                                     const std::string& val) {
  const int len = 31*9;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, xz,
           yx, yy, yz,
           zx, zy, zz;
    sscanf(valptr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf ",
           &xx, &xy, &xz,
           &yx, &yy, &yz,
           &zx, &zy, &zz);
    field(i).xx(xx);
    field(i).xy(xy);
    field(i).xz(xz);
    field(i).yx(yx);
    field(i).yy(yy);
    field(i).yz(yz);
    field(i).zx(zx);
    field(i).zy(zy);
    field(i).zz(zz);
  }
  delete [] valptr;
}

//------------------------------------------------------------------------------
// Specializations: string -> Field<D, SymTensor>.
//------------------------------------------------------------------------------
template<>
inline
void
string2field<Dim<1>, Dim<1>::SymTensor>(FieldSpace::Field<Dim<1>, Dim<1>::SymTensor>& field,
                                        const std::string& val) {
  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx;
    sscanf(valptr, "%lf", &xx);
    field(i).xx(xx);
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<2>, Dim<2>::SymTensor>(FieldSpace::Field<Dim<2>, Dim<2>::SymTensor>& field,
                                        const std::string& val) {
  const int len = 31*3;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, yy;
    sscanf(valptr, "%lf %lf %lf", &xx, &xy, &yy);
    field(i).xx(xx);
    field(i).xy(xy);
    field(i).yy(yy);
  }
  delete [] valptr;
}

template<>
inline
void
string2field<Dim<3>, Dim<3>::SymTensor>(FieldSpace::Field<Dim<3>, Dim<3>::SymTensor>& field,
                                        const std::string& val) {
  const int len = 31*6;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int n;
  sscanf(valptr, "%i", &n);
  VERIFY(field.nodeList().numInternalNodes() == n);

  for (int i = 0; i != n; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, xz,
               yy, yz,
                   zz;
    sscanf(valptr, "%lf %lf %lf %lf %lf %lf ",
           &xx, &xy, &xz,
                &yy, &yz,
                     &zz);
    field(i).xx(xx);
    field(i).xy(xy);
    field(i).xz(xz);
    field(i).yy(yy);
    field(i).yz(yz);
    field(i).zz(zz);
  }
  delete [] valptr;
}

}
