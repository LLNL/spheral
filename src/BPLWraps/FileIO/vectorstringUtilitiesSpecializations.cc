#include <iostream>
#include <iomanip>
#include <sstream>

#include "vectorstringUtilities.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// vector<int> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<int>& val,
              const int precision) {
  std::string result;

  char tmp[21];
  tmp[20] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<int>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%20i ", *itr);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<bool> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<bool>& val,
              const int precision) {
  std::string result;

  char tmp[21];
  tmp[20] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<bool>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%i ", *itr);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<double> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<double>& val,
              const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<double>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e ", *itr);
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<std::string> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<std::string>& val,
              const int precision) {
  std::string result;

  char tmp[21];
  tmp[20] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<std::string>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    char tmpsize[7];
    tmpsize[6] = 0;
    sprintf(tmpsize, "%6i", itr->size());
    result += std::string(tmpsize);

    char tmp[val.size() + 2];
    tmp[val.size() + 1] = 0;
    sprintf(tmp, "%s ", (*itr).c_str());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<Vector1d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<1>::Vector>& val,
              const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<1>::Vector>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e ", itr->x());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<Vector2d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<2>::Vector>& val,
              const int precision) {
  std::string result;

  char tmp[31*2 + 1];
  tmp[31*2] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<2>::Vector>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e %30.20e ", itr->x(), itr->y());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<Vector3d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<3>::Vector>& val,
              const int precision) {
  std::string result;

  char tmp[31*3 + 1];
  tmp[31*3] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<3>::Vector>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e %30.20e %30.20e ", itr->x(), itr->y(), itr->z());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<Tensor1d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<1>::Tensor>& val,
              const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<1>::Tensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e ", itr->xx());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<Tensor2d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<2>::Tensor>& val,
              const int precision) {
  std::string result;

  char tmp[31*4 + 1];
  tmp[31*4] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<2>::Tensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e %30.20e %30.20e %30.20e ",
            itr->xx(), itr->xy(),
            itr->yx(), itr->yy());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<Tensor3d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<3>::Tensor>& val,
              const int precision) {
  std::string result;

  char tmp[31*9 + 1];
  tmp[31*9] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<3>::Tensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e %30.20e ",
            itr->xx(), itr->xy(), itr->xz(),
            itr->yx(), itr->yy(), itr->yz(),
            itr->zx(), itr->zy(), itr->zz());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<SymTensor1d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<1>::SymTensor>& val,
              const int precision) {
  std::string result;

  char tmp[32];
  tmp[31] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<1>::SymTensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e ", itr->xx());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<SymTensor2d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<2>::SymTensor>& val,
              const int precision) {
  std::string result;

  char tmp[31*3 + 1];
  tmp[31*3] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<2>::SymTensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e %30.20e %30.20e ",
            itr->xx(), itr->xy(),
                       itr->yy());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<SymTensor3d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<3>::SymTensor>& val,
              const int precision) {
  std::string result;

  char tmp[31*6 + 1];
  tmp[31*6] = 0;
  sprintf(tmp, "%20i", val.size());
  result += std::string(tmp);

  for (std::vector<Spheral::Dim<3>::SymTensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) {
    sprintf(tmp, "%30.20e %30.20e %30.20e %30.20e %30.20e %30.20e ",
            itr->xx(), itr->xy(), itr->xz(),
                       itr->yy(), itr->yz(),
                                  itr->zz());
    result += std::string(tmp);
  }
  result += " ";

  return result;
}

//------------------------------------------------------------------------------
// vector<ThirdRankTensor1d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<1>::ThirdRankTensor>& val,
              const int precision) {
  typedef Spheral::Dim<1>::ThirdRankTensor Element;
  std::stringstream result;

  // Output the size of our vector.
  // result << std::setw(30);
  result << val.size() << " ";

  // Now all the individual elements.
  result << std::setprecision(20) << std::scientific;
  for (std::vector<Spheral::Dim<1>::ThirdRankTensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) result << *itr << " ";

  return result.str();
}

//------------------------------------------------------------------------------
// vector<ThirdRankTensor2d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<2>::ThirdRankTensor>& val,
              const int precision) {
  typedef Spheral::Dim<2>::ThirdRankTensor Element;
  std::stringstream result;

  // Output the size of our vector.
  // result << std::setw(30);
  result << val.size() << " ";

  // Now all the individual elements.
  result << std::setprecision(20) << std::scientific;
  for (std::vector<Spheral::Dim<2>::ThirdRankTensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) result << *itr << " ";

  return result.str();
}

//------------------------------------------------------------------------------
// vector<ThirdRankTensor3d> -> string conversion.
//------------------------------------------------------------------------------
template<>
std::string
vector2string(const std::vector<Spheral::Dim<3>::ThirdRankTensor>& val,
              const int precision) {
  typedef Spheral::Dim<3>::ThirdRankTensor Element;
  std::stringstream result;

  // Output the size of our vector.
  // result << std::setw(30);
  result << val.size() << " ";

  // Now all the individual elements.
  result << std::setprecision(20) << std::scientific;
  for (std::vector<Spheral::Dim<3>::ThirdRankTensor>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) result << *itr << " ";

  return result.str();
}

//------------------------------------------------------------------------------
// string -> vector<int> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<int>
string2vector<int>(const std::string& val) {
  std::vector<int> result;

  const int len = 21;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, 20, 20 + len*i);
    sscanf(valptr, "%i ", &result[i]);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<bool> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<bool>
string2vector<bool>(const std::string& val) {
  std::vector<bool> result;

  const int len = 2;
  char* valptr = new char[21];
  valptr[20] = 0;

  val.copy(valptr, 20, 0);
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  valptr[1] = 0;
  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    int tmp;
    sscanf(valptr, "%i", &tmp);
    result[i] = tmp == 1 ? true : false;
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<double> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<double>
string2vector<double>(const std::string& val) {
  std::vector<double> result;

  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    sscanf(valptr, "%lf", &result[i]);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<string> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<std::string>
string2vector<std::string>(const std::string& val) {
  std::vector<std::string> result;

  char* valptr = new char[21];
  valptr[20] = 0;

  val.copy(valptr, 20, 0);
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);
  delete [] valptr;

  int offset = 20;
  for (int i = 0; i != size; ++i) {
    VERIFY(offset + 6 <= val.size());
    valptr = new char[7];
    valptr[6] = 0;
    val.copy(valptr, 6, offset);
    int size;
    sscanf(valptr, "%i", &size);
    delete [] valptr;

    VERIFY(size >= 0);
    VERIFY(offset + 6 + size + 1 <= val.size());
    valptr = new char[size + 1];
    valptr[size] = 0;
    val.copy(valptr, size, offset + 6);
    result[i] = std::string(valptr);
    offset += 6 + size + 1;
    delete [] valptr;
  }

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<Vector1d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<1>::Vector>
string2vector<Spheral::Dim<1>::Vector>(const std::string& val) {
  std::vector<Spheral::Dim<1>::Vector> result;

  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double x;
    sscanf(valptr, "%lf", &x);
    result[i].x(x);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<Vector2d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<2>::Vector>
string2vector<Spheral::Dim<2>::Vector>(const std::string& val) {
  std::vector<Spheral::Dim<2>::Vector> result;

  const int len = 31*2;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double x, y;
    sscanf(valptr, "%lf %lf", &x, &y);
    result[i].x(x);
    result[i].y(y);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<Vector3d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<3>::Vector>
string2vector<Spheral::Dim<3>::Vector>(const std::string& val) {
  std::vector<Spheral::Dim<3>::Vector> result;

  const int len = 31*3;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double x, y, z;
    sscanf(valptr, "%lf %lf %lf", &x, &y, &z);
    result[i].x(x);
    result[i].y(y);
    result[i].z(z);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<Tensor1d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<1>::Tensor>
string2vector<Spheral::Dim<1>::Tensor>(const std::string& val) {
  std::vector<Spheral::Dim<1>::Tensor> result;

  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx;
    sscanf(valptr, "%lf", &xx);
    result[i].xx(xx);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<Tensor2d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<2>::Tensor>
string2vector<Spheral::Dim<2>::Tensor>(const std::string& val) {
  std::vector<Spheral::Dim<2>::Tensor> result;

  const int len = 31*4;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, 
           yx, yy;
    sscanf(valptr, "%lf %lf %lf %lf", &xx, &xy, &yx, &yy);
    result[i].xx(xx);
    result[i].xy(xy);
    result[i].yx(yx);
    result[i].yy(yy);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<Tensor3d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<3>::Tensor>
string2vector<Spheral::Dim<3>::Tensor>(const std::string& val) {
  std::vector<Spheral::Dim<3>::Tensor> result;

  const int len = 31*9;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, xz,
           yx, yy, yz,
           zx, zy, zz;
    sscanf(valptr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &xx, &xy, &xz,
           &yx, &yy, &yz,
           &zx, &zy, &zz);
    result[i].xx(xx);
    result[i].xy(xy);
    result[i].xz(xz);
    result[i].yx(yx);
    result[i].yy(yy);
    result[i].yz(yz);
    result[i].zx(zx);
    result[i].zy(zy);
    result[i].zz(zz);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<SymTensor1d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<1>::SymTensor>
string2vector<Spheral::Dim<1>::SymTensor>(const std::string& val) {
  std::vector<Spheral::Dim<1>::SymTensor> result;

  const int len = 31;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx;
    sscanf(valptr, "%lf", &xx);
    result[i].xx(xx);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<SymTensor2d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<2>::SymTensor>
string2vector<Spheral::Dim<2>::SymTensor>(const std::string& val) {
  std::vector<Spheral::Dim<2>::SymTensor> result;

  const int len = 31*3;
  char* valptr = new char[len + 1];
  valptr[len] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, 
               yy;
    sscanf(valptr, "%lf %lf %lf", &xx, &xy, &yy);
    result[i].xx(xx);
    result[i].xy(xy);
    result[i].yy(yy);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<SymTensor3d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<3>::SymTensor>
string2vector<Spheral::Dim<3>::SymTensor>(const std::string& val) {
  std::vector<Spheral::Dim<3>::SymTensor> result;

  const int len = 31*6;
  char* valptr = new char[len + 1];
  valptr[31*6] = 0;

  val.copy(valptr, 20, 0);
  valptr[20] = 0;
  int size;
  sscanf(valptr, "%i", &size);
  result.resize(size);
  VERIFY(result.size() == size);

  for (int i = 0; i != size; ++i) {
    VERIFY(20 + len*(i + 1) <= val.size());
    val.copy(valptr, len, 20 + len*i);
    double xx, xy, xz,
               yy, yz,
                   zz;
    sscanf(valptr, "%lf %lf %lf %lf %lf %lf ",
           &xx, &xy, &xz,
                &yy, &yz,
                     &zz);
    result[i].xx(xx);
    result[i].xy(xy);
    result[i].xz(xz);
    result[i].yy(yy);
    result[i].yz(yz);
    result[i].zz(zz);
  }
  delete [] valptr;

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<ThirdRankTensor1d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<1>::ThirdRankTensor>
string2vector<Spheral::Dim<1>::ThirdRankTensor>(const std::string& val) {

  typedef Spheral::Dim<1>::ThirdRankTensor Element;
  std::vector<Spheral::Dim<1>::ThirdRankTensor> result;

  // Prepare to read from the string.
  std::istringstream is(val);

  // Get the size of the vector.
  int size;
  is >> size;
  result.reserve(size);

  // Now read the individual values.
  for (size_t i = 0; i != size; ++i) {
    Element val;
    is >> val;
    result.push_back(val);
  }

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<ThirdRankTensor2d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<2>::ThirdRankTensor>
string2vector<Spheral::Dim<2>::ThirdRankTensor>(const std::string& val) {

  typedef Spheral::Dim<2>::ThirdRankTensor Element;
  std::vector<Spheral::Dim<2>::ThirdRankTensor> result;

  // Prepare to read from the string.
  std::istringstream is(val);

  // Get the size of the vector.
  int size;
  is >> size;
  result.reserve(size);

  // Now read the individual values.
  for (size_t i = 0; i != size; ++i) {
    Element val;
    is >> val;
    result.push_back(val);
  }

  return result;
}

//------------------------------------------------------------------------------
// string -> vector<ThirdRankTensor3d> conversion.
//------------------------------------------------------------------------------
template<>
std::vector<Spheral::Dim<3>::ThirdRankTensor>
string2vector<Spheral::Dim<3>::ThirdRankTensor>(const std::string& val) {

  typedef Spheral::Dim<3>::ThirdRankTensor Element;
  std::vector<Spheral::Dim<3>::ThirdRankTensor> result;

  // Prepare to read from the string.
  std::istringstream is(val);

  // Get the size of the vector.
  int size;
  is >> size;
  result.reserve(size);

  // Now read the individual values.
  for (size_t i = 0; i != size; ++i) {
    Element val;
    is >> val;
    result.push_back(val);
  }

  return result;
}

}
