#ifndef __PBGWRAPS_SILOWRAPPERS__
#define __PBGWRAPS_SILOWRAPPERS__

#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include "memory"

#include "silo.h"

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace silo {

//------------------------------------------------------------------------------
// Trait class to help mapping Spheral types to silo.
//------------------------------------------------------------------------------
template<typename T>
struct Spheral2Silo;

template<>
struct Spheral2Silo<double> {
  typedef double Value;
  static unsigned numElements() { return 1; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x;
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::Vector> {
  typedef Spheral::Dim<2>::Vector Value;
  static unsigned numElements() { return 2; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.x();
    silovars[1][i] = x.y();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::Tensor> {
  typedef Spheral::Dim<2>::Tensor Value;
  static unsigned numElements() { return 4; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy();
    silovars[2][i] = x.yx(); silovars[3][i] = x.yy();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<2>::SymTensor> {
  typedef Spheral::Dim<2>::SymTensor Value;
  static unsigned numElements() { return 4; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy();
    silovars[2][i] = x.yx(); silovars[3][i] = x.yy();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::Vector> {
  typedef Spheral::Dim<3>::Vector Value;
  static unsigned numElements() { return 3; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.x();
    silovars[1][i] = x.y();
    silovars[2][i] = x.z();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::Tensor> {
  typedef Spheral::Dim<3>::Tensor Value;
  static unsigned numElements() { return 9; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy(); silovars[2][i] = x.xz();
    silovars[3][i] = x.yx(); silovars[4][i] = x.yy(); silovars[5][i] = x.yz();
    silovars[6][i] = x.zx(); silovars[7][i] = x.zy(); silovars[8][i] = x.zz();
  }
};

template<>
struct Spheral2Silo<Spheral::Dim<3>::SymTensor> {
  typedef Spheral::Dim<3>::SymTensor Value;
  static unsigned numElements() { return 9; }
  static void copyElement(const Value& x, double** silovars, const unsigned i) {
    silovars[0][i] = x.xx(); silovars[1][i] = x.xy(); silovars[2][i] = x.xz();
    silovars[3][i] = x.yx(); silovars[4][i] = x.yy(); silovars[5][i] = x.yz();
    silovars[6][i] = x.zx(); silovars[7][i] = x.zy(); silovars[8][i] = x.zz();
  }
};

//------------------------------------------------------------------------------
// Names
//------------------------------------------------------------------------------
// typedef ::DBfile DBfile;
// typedef ::DBoptlist DBoptlist;
// typedef ::DBmrgtree DBmrgtree;
// using ::DBMakeOptlist;
// using ::DBClearOptlist;
// using ::DBClearOption;
// using std::string;

//------------------------------------------------------------------------------
// Convert a std::string -> char*
//------------------------------------------------------------------------------
struct ConvertStringToCharStar {
  char* operator()(const std::string& x) {
    return const_cast<char*>(x.c_str());
  }
};

//------------------------------------------------------------------------------
// A trait class for for mapping types -> silo types.
//------------------------------------------------------------------------------
template<typename T> struct SiloTraits {};

template<>
struct SiloTraits<int> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_INT; }
};

template<>
struct SiloTraits<short> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_SHORT; }
};

template<>
struct SiloTraits<long> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_LONG; }
};

template<>
struct SiloTraits<float> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_FLOAT; }
};

template<>
struct SiloTraits<double> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_DOUBLE; }
};

template<>
struct SiloTraits<char> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_CHAR; }
};

template<>
struct SiloTraits<long long> {
  static std::vector<int> dims() { return std::vector<int>(1, 1); }
  static int datatype() { return DB_LONG_LONG; }
};

//------------------------------------------------------------------------------
// Wrapper class to handle the memory managemnt necessary with DBoptlist.
//------------------------------------------------------------------------------
struct DBoptlist_wrapper {
  DBoptlist* mOptlistPtr;
  std::vector<std::shared_ptr<void> > mCache;

  // Constructors.
  DBoptlist_wrapper(const int maxopts=1024):
    mOptlistPtr(DBMakeOptlist(maxopts)),
    mCache() {}

  // Destructor.
  ~DBoptlist_wrapper() {
    VERIFY(DBFreeOptlist(mOptlistPtr) == 0);
  }

  // Generic functor definitions for adding and getting options.
  template<typename Value>
  struct AddOptionFunctor {
    int writeValue(DBoptlist_wrapper& optlist_wrapper,
                   const int option,
                   const Value& value) {
      std::shared_ptr<void> voidValue(new Value(value));
      optlist_wrapper.mCache.push_back(voidValue);
      return DBAddOption(optlist_wrapper.mOptlistPtr, option, voidValue.get());
    }
    int writeVector(DBoptlist_wrapper& optlist_wrapper,
                    const int option,
                    const int option_size,
                    const std::vector<Value>& value) {
      DBoptlist_wrapper::AddOptionFunctor<int>().writeValue(optlist_wrapper, option_size, value.size());
      std::shared_ptr<void> voidValue(new std::vector<Value>(value));
      optlist_wrapper.mCache.push_back(voidValue);
      Value* frontPtr = &(((std::vector<Value>*) voidValue.get())->front());
      return DBAddOption(optlist_wrapper.mOptlistPtr, option, frontPtr);
    }
  };

  template<typename Value>
  struct GetOptionFunctor {
    Value readValue(DBoptlist_wrapper& optlist_wrapper,
                    const int option) {
      return *((Value*) DBGetOption(optlist_wrapper.mOptlistPtr, option));
    }
    std::vector<Value> readVector(DBoptlist_wrapper& optlist_wrapper,
                                  const int option,
                                  const int option_size) {
      const unsigned vecsize = DBoptlist_wrapper::GetOptionFunctor<int>().readValue(optlist_wrapper, option_size);
      Value* frontPtr = (Value*) DBGetOption(optlist_wrapper.mOptlistPtr, option);
      return std::vector<Value>(frontPtr, frontPtr + vecsize);
    }
  };

  // Function definitions that use the functors.
  template<typename Value>
  int addOption(const int option,
                const Value& value) {
    return DBoptlist_wrapper::AddOptionFunctor<Value>().writeValue(*this, option, value);
  }

  template<typename Value>
  int addOption(const int option,
                const int option_size,
                const std::vector<Value>& value) {
    return DBoptlist_wrapper::AddOptionFunctor<Value>().writeVector(*this, option, option_size, value);
  }

  template<typename Value>
  Value getOption(const int option) {
    return DBoptlist_wrapper::GetOptionFunctor<Value>().readValue(*this, option);
  }

  template<typename Value>
  std::vector<Value> getOption(const int option,
                               const int option_size) {
    return DBoptlist_wrapper::GetOptionFunctor<Value>().readVector(*this, option, option_size);
  }
};
                              
//------------------------------------------------------------------------------
// Wrapper class to handle the memory managemnt necessary with DBmrgtree
//------------------------------------------------------------------------------
struct DBmrgtree_wrapper {
  DBmrgtree* mDBmrgtree;

  // Constructors.
  DBmrgtree_wrapper(int mesh_type,
                    int info_bits,
                    int max_children,
                    DBoptlist_wrapper optlist):
    mDBmrgtree(DBMakeMrgtree(mesh_type, info_bits, max_children, optlist.mOptlistPtr)) {}

  // Destructor.
  ~DBmrgtree_wrapper() {
    DBFreeMrgtree(mDBmrgtree);
  }

  // name
  std::string name() const { return (mDBmrgtree->name != NULL ? 
                                     std::string(mDBmrgtree->name) :
                                     std::string()); }
  void name(std::string val) { 
    mDBmrgtree->name = new char[val.length() + 1];
    strcpy(mDBmrgtree->name, val.c_str());
  }

  // src_mesh_name
  std::string src_mesh_name() const { return (mDBmrgtree->src_mesh_name != NULL ? 
                                              std::string(mDBmrgtree->src_mesh_name) :
                                              std::string()); }
  void src_mesh_name(std::string val) { 
    mDBmrgtree->src_mesh_name = new char[val.length() + 1];
    strcpy(mDBmrgtree->src_mesh_name, val.c_str());
  }

  // src_mesh_type
  int src_mesh_type() const { return mDBmrgtree->src_mesh_type; }
  void src_mesh_type(int val) { mDBmrgtree->src_mesh_type = val; }

  // type_info_bits
  int type_info_bits() const { return mDBmrgtree->type_info_bits; }
  void type_info_bits(int val) { mDBmrgtree->type_info_bits = val; }

  // num_nodes
  int num_nodes() const { return mDBmrgtree->num_nodes; }
  void num_nodes(int val) { mDBmrgtree->num_nodes = val; }
};
                              
//..............................................................................
// std::string specializations.
//..............................................................................
template<>
struct
DBoptlist_wrapper::AddOptionFunctor<std::string> {
  int
  writeValue(DBoptlist_wrapper& optlist_wrapper,
             const int option,
             const std::string& value) {
    std::shared_ptr<void> voidValue(new std::string(value));
    optlist_wrapper.mCache.push_back(voidValue);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, (char*) ((std::string*) voidValue.get())->c_str());
  }
  int
  writeVector(DBoptlist_wrapper& optlist_wrapper,
              const int option,
              const int option_size,
              const std::vector<std::string>& value) {
    VERIFY(optlist_wrapper.addOption<int>(option_size, value.size()) == 0);
    std::shared_ptr<void> voidCopy(new std::vector<std::string>(value));
    std::shared_ptr<void> voidValue(new std::vector<char*>());
    std::vector<std::string>& stringvec = *((std::vector<std::string>*) voidCopy.get());
    std::vector<char*>& charvec = *((std::vector<char*>*) (voidValue.get()));
    for (unsigned k = 0; k != value.size(); ++k) charvec.push_back(const_cast<char*>(stringvec[k].c_str()));
    VERIFY(charvec.size() == value.size());
    optlist_wrapper.mCache.push_back(voidValue);
    return DBAddOption(optlist_wrapper.mOptlistPtr, option, (char**) &charvec.front());
  }
};

template<>
struct
DBoptlist_wrapper::GetOptionFunctor<std::string> {
  std::string
  readValue(DBoptlist_wrapper& optlist_wrapper,
            const int option) {
    char* result = (char*) DBGetOption(optlist_wrapper.mOptlistPtr, option);
    return std::string(result);
  }
  std::vector<std::string>
  readVector(DBoptlist_wrapper& optlist_wrapper,
             const int option,
             const int option_size) {
    const int resultsize = optlist_wrapper.getOption<int>(option_size);
    VERIFY(resultsize > 0);
    char** chararray = (char**) DBGetOption(optlist_wrapper.mOptlistPtr, option);
    std::vector<std::string> result(chararray, chararray + resultsize);
    // for (unsigned k = 0;  k != resultsize; ++k) result.push_back(std::string(chararray[k]));
    VERIFY(result.size() == resultsize);
    return result;
  }
};

//------------------------------------------------------------------------------
// DBCreate
//------------------------------------------------------------------------------
inline
DBfile*
DBCreate_wrap(std::string pathName,
              int mode,
              int target,
              std::string fileInfo,
              int fileType) {
  return DBCreate(pathName.c_str(), mode, target, fileInfo.c_str(), fileType);
}

//------------------------------------------------------------------------------
// DBOpen
//------------------------------------------------------------------------------
inline
DBfile*
DBOpen_wrap(std::string pathName,
              int type,
              int mode) {
  return DBOpen(pathName.c_str(), type, mode);
}

//------------------------------------------------------------------------------
// DBMakeMrgtree
//------------------------------------------------------------------------------
// inline
// DBmrgtree*
// DBMakeMrgtree_wrap(int mesh_type,
//                    int info_bits,
//                    int max_children,
//                    DBoptlist_wrapper& optlist) {
//   return DBMakeMrgtree(mesh_type,
//                        info_bits,
//                        max_children,
//                        optlist.mOptlistPtr);
// }

//------------------------------------------------------------------------------
// DBFreeMrgtree
//------------------------------------------------------------------------------
// inline
// void
// DBFreeMrgtree_wrap(DBmrgtree& tree) {
//   DBFreeMrgtree(&tree);
// }

//------------------------------------------------------------------------------
// DBClose
//------------------------------------------------------------------------------
inline
int
DBClose(DBfile& file) {
  return DBClose(&file);
}

//------------------------------------------------------------------------------
// DBWrite
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBWrite(DBfile& file,
        std::string varname,
        T& var) {
  return DBWrite(&file, varname.c_str(), (void*) &var, 
                 &(SiloTraits<T>::dims()).front(),
                 SiloTraits<T>::dims().size(),
                 SiloTraits<T>::datatype());
}

//------------------------------------------------------------------------------
// DBReadVar
//------------------------------------------------------------------------------
template<typename T>
inline
T
DBReadVar(DBfile& file,
          std::string varname) {
  T result;
  DBReadVar(&file, varname.c_str(), (void*) &result);
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultimesh
//------------------------------------------------------------------------------
inline
int
DBPutMultimesh(DBfile& file,
               std::string name,
               std::vector<std::string>& meshNames,
               std::vector<int>& meshTypes,
               DBoptlist_wrapper& optlist) {

  // Pre-conditions.
  VERIFY2(meshNames.size() == meshTypes.size(), "meshNames and meshTypes must be same length:  " << meshNames.size() << " != " << meshTypes.size());

  // Convert names to char*.
  std::vector<char*> meshNames1;
  std::transform(meshNames.begin(), meshNames.end(), std::back_inserter(meshNames1),
                 ConvertStringToCharStar());
  CHECK(meshNames1.size() == meshNames.size());

  // Do the deed.
  const int result = DBPutMultimesh(&file, 
                                    name.c_str(), 
                                    meshNames.size(),
                                    &meshNames1.front(), 
                                    &meshTypes.front(),
                                    optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultimat
//------------------------------------------------------------------------------
inline
int
DBPutMultimat(DBfile& file,
              std::string name,
              std::vector<std::string>& matNames,
              DBoptlist_wrapper& optlist) {

  // Convert names to char*.
  std::vector<char*> matNames1;
  std::transform(matNames.begin(), matNames.end(), std::back_inserter(matNames1),
                 ConvertStringToCharStar());
  CHECK(matNames1.size() == matNames.size());

  // Do the deed.
  const int result = DBPutMultimat(&file, 
                                   name.c_str(), 
                                   matNames.size(),
                                   &matNames1.front(), 
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutCompoundarray
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutCompoundarray(DBfile& file,
                   std::string name,
                   std::vector<std::string>& elemNames,
                   std::vector<std::vector<T> >& values,
                   DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(elemNames.size() == values.size());

  // Convert names to char*.
  std::vector<char*> elemNames1;
  std::transform(elemNames.begin(), elemNames.end(), std::back_inserter(elemNames1),
                 ConvertStringToCharStar());
  CHECK(elemNames1.size() == elemNames.size());

  // Read the sizes of each array.
  std::vector<int> elemLengths;
  elemLengths.reserve(values.size());
  for (unsigned k = 0; k != values.size(); ++k) elemLengths.push_back(values[k].size());
  const unsigned numValues = std::accumulate(elemLengths.begin(), elemLengths.end(), 0);

  // Flatten the values to a single arrray.
  std::vector<T> flatValues;
  flatValues.reserve(numValues);
  for (unsigned k = 0; k != values.size(); ++k) {
    for (unsigned j = 0; j != values[k].size(); ++j) {
      flatValues.push_back(values[k][j]);
    }
  }
  CHECK(flatValues.size() == numValues);

  // Do the deed.
  CHECK(elemNames.size() == elemLengths.size());
  const int result = DBPutCompoundarray(&file, 
                                        name.c_str(), 
                                        &elemNames1.front(),
                                        &elemLengths.front(),
                                        elemNames1.size(),
                                        &flatValues.front(),
                                        numValues,
                                        SiloTraits<T>::datatype(),
                                        optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMultivar
//------------------------------------------------------------------------------
inline
int
DBPutMultivar(DBfile& file,
              std::string name,
              std::vector<std::string>& varNames,
              std::vector<int>& varTypes,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned numVars = varNames.size();
  VERIFY(varTypes.size() == numVars);

  // Convert names to char*.
  std::vector<char*> varNames1;
  std::transform(varNames.begin(), varNames.end(), std::back_inserter(varNames1),
                 ConvertStringToCharStar());
  CHECK(varNames1.size() == varNames.size());

  // Do the deed.
  const int result = DBPutMultivar(&file, 
                                   const_cast<char*>(name.c_str()), 
                                   numVars,
                                   &varNames1.front(),
                                   &varTypes.front(),
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutMaterial
//------------------------------------------------------------------------------
inline
int
DBPutMaterial(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<int>& matnos,
              std::vector<int>& matlist,
              std::vector<int>& mix_next,
              std::vector<int>& mix_mat,
              std::vector<int>& mix_zone,
              std::vector<double>& mix_vf,
              DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nmat = matnos.size();
  const unsigned numMix = mix_next.size();
  VERIFY(mix_mat.size() == numMix);
  VERIFY(mix_zone.size() == numMix);
  VERIFY(mix_vf.size() == numMix);

  // Dimensionality of the matlist.
  // For now we only support 1-D lists (I don't understand what silo does with 
  // greater dimensionality?)
  std::vector<int> dims(1, matlist.size());

  // Do the deed.
  const int result = DBPutMaterial(&file, 
                                   name.c_str(),
                                   meshName.c_str(),
                                   nmat,
                                   &matnos.front(),
                                   &matlist.front(),
                                   &dims.front(),
                                   dims.size(),
                                   &mix_next.front(),
                                   &mix_mat.front(),
                                   &mix_zone.front(),
                                   &mix_vf.front(),
                                   numMix,
                                   SiloTraits<double>::datatype(),
                                   optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutUcdmesh
//------------------------------------------------------------------------------
inline
int
DBPutUcdmesh(DBfile& file,
             std::string name,
             std::vector<std::vector<double> >& coords,
             int nzones,
             std::string zonel_name,
             std::string facel_name,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);
  const unsigned nnodes = coords[0].size();
  for (unsigned idim = 0; idim != ndims; ++idim) VERIFY(coords[idim].size() == nnodes);
  VERIFY(nzones > 0);

  // We need the C-stylish pointers to the coordinates.
  double** coordPtrs = new double*[ndims];
  for (unsigned k = 0; k != ndims; ++k) {
    coordPtrs[k] = new double[nnodes];
    std::copy(coords[k].begin(), coords[k].end(), coordPtrs[k]);
  }

  // Convert strings to char*.
  char* zonel_name1 = (zonel_name == "NULL") ? NULL : const_cast<char*>(zonel_name.c_str());
  char* facel_name1 = (facel_name == "NULL") ? NULL : const_cast<char*>(facel_name.c_str());

  // Do the deed.
  const int result = DBPutUcdmesh(&file, 
                                  name.c_str(),
                                  ndims,
                                  NULL,
                                  coordPtrs,
                                  nnodes,
                                  nzones,
                                  zonel_name1,
                                  facel_name1,
                                  SiloTraits<double>::datatype(),
                                  optlist.mOptlistPtr);

  // That's it.
  for (unsigned k = 0; k != ndims; ++k) delete[] coordPtrs[k];
  delete[] coordPtrs;
  return result;
}

//------------------------------------------------------------------------------
// DBPutDefvars
//------------------------------------------------------------------------------
inline
int
DBPutDefvars(DBfile& file,
             std::string name,
             std::vector<std::string>& varNames,
             std::vector<int>& varTypes,
             std::vector<std::string>& varDefs,
             std::vector<DBoptlist_wrapper*>& optlists) {

  // Preconditions.
  const unsigned ndefs = varNames.size();
  VERIFY(varNames.size() == ndefs and
         varTypes.size() == ndefs and
         varDefs.size() == ndefs and
         optlists.size() == ndefs);

  // Convert names to char*'s
  std::vector<char*> names, defns;
  std::transform(varNames.begin(), varNames.end(), std::back_inserter(names), ConvertStringToCharStar());
  std::transform(varDefs.begin(), varDefs.end(), std::back_inserter(defns), ConvertStringToCharStar());
  CHECK(names.size() == varNames.size());
  CHECK(defns.size() == varDefs.size());
  
  // Copy the optlists to an array of pointers.
  std::vector<DBoptlist*> optlistptrs;
  for (std::vector<DBoptlist_wrapper*>::iterator itr = optlists.begin();
       itr != optlists.end();
       ++itr) optlistptrs.push_back((*itr)->mOptlistPtr);
  CHECK(optlistptrs.size() == ndefs);

  return DBPutDefvars(&file,
                      name.c_str(),
                      ndefs,
                      &names.front(),
                      &varTypes.front(),
                      &defns.front(),
                      &optlistptrs.front());
}

//------------------------------------------------------------------------------
// DBPutUcdvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutUcdvar(DBfile& file,
            std::string name,
            std::string meshName,
            std::vector<T>& values,
            std::vector<T>& mixValues,
            int centering,
            DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);

  unsigned i, j;

  // Build the sub-variable names.
  int nvars = Spheral2Silo<T>::numElements();
  int nels = values.size();
  int mixlen = mixValues.size();
  std::vector<char*> varnames;
  for (i = 0; i != nvars; ++i) {
    varnames.push_back(const_cast<char*>((name + "_").c_str()));
    sprintf(varnames.back(), "%i", i);
  }

  // Build the sub-variables.
  double** vars = new double*[nvars];
  double** mixvars = new double*[nvars];
  for (i = 0; i != nvars; ++i) {
    vars[i] = new double[nels];
    mixvars[i] = new double[mixlen];
  }
  for (j = 0; j != nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);
  for (j = 0; j != mixlen; ++j) Spheral2Silo<T>::copyElement(mixValues[j], mixvars, j);

  const int result = DBPutUcdvar(&file,
                                 name.c_str(),
                                 meshName.c_str(),
                                 nvars,
                                 &varnames.front(),
                                 (void*) vars,
                                 nels,
                                 (void*) mixvars,
                                 mixlen,
                                 SiloTraits<double>::datatype(),
                                 centering,
                                 optlist.mOptlistPtr);

  // That's it.
  for (i = 0; i != nvars; ++i) {
    delete[] vars[i];
    delete[] mixvars[i];
  }
  delete[] vars;
  delete[] mixvars;
  return result;
}

//------------------------------------------------------------------------------
// DBPutUcdvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutUcdvar1(DBfile& file,
             std::string name,
             std::string meshName,
             std::vector<T>& values,
             std::vector<T>& mixValues,
             int centering,
             DBoptlist_wrapper& optlist) {

  // Preconditions.
  VERIFY(centering == DB_NODECENT or
         centering == DB_EDGECENT or
         centering == DB_FACECENT or
         centering == DB_ZONECENT);

  return DBPutUcdvar1(&file,
                      name.c_str(),
                      meshName.c_str(),
                      (void*) &(*values.begin()),
                      values.size(),
                      (void*) &(*mixValues.begin()),
                      mixValues.size(),
                      SiloTraits<T>::datatype(),
                      centering,
                      optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBPutZonelist2
//------------------------------------------------------------------------------
inline
int
DBPutZonelist2(DBfile& file,
               std::string name,
               unsigned ndims,
               std::vector<std::vector<int> >& zoneNodes,
               unsigned low_offset,
               unsigned high_offset,
               std::vector<int>& shapetype,
               std::vector<int>& shapesize,
               std::vector<int>& shapecount,
               DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nzones = zoneNodes.size();
  const unsigned nshapes = shapetype.size();
  VERIFY(shapetype.size() <= nzones);
  VERIFY(shapesize.size() <= nzones);
  VERIFY(shapecount.size() <= nzones);
  VERIFY(shapesize.size() == nshapes);
  VERIFY(shapecount.size() == nshapes);

  // Construct the flat array of zone nodes.
  std::vector<int> nodelist;
  for (unsigned k = 0; k != zoneNodes.size(); ++k) {
    //nodelist.push_back(zoneNodes[k].size());
    std::copy(zoneNodes[k].begin(), zoneNodes[k].end(), std::back_inserter(nodelist));
  }

  // Do the deed.
  const int result = DBPutZonelist2(&file, 
                                    name.c_str(),
                                    nzones,
                                    ndims,
                                    &nodelist.front(),
                                    nodelist.size(),
                                    0,
                                    low_offset,
                                    high_offset,
                                    &shapetype.front(),
                                    &shapesize.front(),
                                    &shapecount.front(),
                                    nshapes,
                                    optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutPHZonelist
//------------------------------------------------------------------------------
inline
int
DBPutPHZonelist(DBfile& file,
                std::string name,
                std::vector<std::vector<int> >& faceNodeLists,
                std::vector<std::vector<int> >& zoneFaceLists,
                unsigned low_offset,
                unsigned high_offset,
                DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned nfaces = faceNodeLists.size();
  const unsigned nzones = zoneFaceLists.size();

  // Construct the flat arrays of face-node info and zone-face info.
  std::vector<int> nodecnts, nodelist, facecnts, facelist;
  for (unsigned k = 0; k != nfaces; ++k) {
    const std::vector<int>& faceNodes = faceNodeLists[k];
    nodecnts.push_back(faceNodes.size());
    std::copy(faceNodes.begin(), faceNodes.end(), std::back_inserter(nodelist));
  }
  for (unsigned k = 0; k != nzones; ++k) {
    const std::vector<int>& zoneFaces = zoneFaceLists[k];
    facecnts.push_back(zoneFaces.size());
    std::copy(zoneFaces.begin(), zoneFaces.end(), std::back_inserter(facelist));
  }
  CHECK(nodecnts.size() == nfaces);
  CHECK(facecnts.size() == nzones);

  // Do the deed.
  const int result = DBPutPHZonelist(&file, 
                                     name.c_str(),
                                     nfaces,
                                     &nodecnts.front(),
                                     nodelist.size(),
                                     &nodelist.front(),
                                     NULL,
                                     nzones,
                                     &facecnts.front(),
                                     facelist.size(),
                                     &facelist.front(),
                                     0,
                                     low_offset,
                                     high_offset,
                                     optlist.mOptlistPtr);

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointmesh
//------------------------------------------------------------------------------
inline
int
DBPutPointmesh(DBfile& file,
               std::string name,
               std::vector<std::vector<double> >& coords,
               DBoptlist_wrapper& optlist) {

  // Preconditions.
  const unsigned ndims = coords.size();
  VERIFY(ndims == 2 or ndims == 3);
  const unsigned npoints = coords[0].size();
  for (unsigned idim = 0; idim != ndims; ++idim) VERIFY(coords[idim].size() == npoints);

  // We need the C-stylish pointers to the coordinates.
  double** coordPtrs = new double*[ndims];
  for (unsigned k = 0; k != ndims; ++k) {
    coordPtrs[k] = new double[npoints];
    std::copy(coords[k].begin(), coords[k].end(), coordPtrs[k]);
  }

  // Do the deed.
  const int result = DBPutPointmesh(&file, 
                                    name.c_str(),
                                    ndims,
                                    coordPtrs,
                                    npoints,
                                    SiloTraits<double>::datatype(),
                                    optlist.mOptlistPtr);

  // That's it.
  for (unsigned k = 0; k != ndims; ++k) delete[] coordPtrs[k];
  delete[] coordPtrs;
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointvar
// We assume here that the underlying element type is double.
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutPointvar(DBfile& file,
              std::string name,
              std::string meshName,
              std::vector<T>& values,
              DBoptlist_wrapper& optlist) {

  unsigned i, j;

  // Build the sub-variable names.
  int nvars = Spheral2Silo<T>::numElements();
  int nels = values.size();

  // Build the sub-variables.
  double** vars = new double*[nvars];
  for (i = 0; i != nvars; ++i) {
    vars[i] = new double[nels];
  }
  for (j = 0; j != nels; ++j) Spheral2Silo<T>::copyElement(values[j], vars, j);

  const int result = DBPutPointvar(&file,
                                   name.c_str(),
                                   meshName.c_str(),
                                   nvars,
                                   (void*) vars,
                                   nels,
                                   SiloTraits<double>::datatype(),
                                   optlist.mOptlistPtr);

  // That's it.
  for (i = 0; i != nvars; ++i) {
    delete[] vars[i];
  }
  delete[] vars;
  return result;
}

//------------------------------------------------------------------------------
// DBPutPointvar1
//------------------------------------------------------------------------------
template<typename T>
inline
int
DBPutPointvar1(DBfile& file,
               std::string name,
               std::string meshName,
               std::vector<T>& values,
               DBoptlist_wrapper& optlist) {

  return DBPutPointvar1(&file,
                        name.c_str(),
                        meshName.c_str(),
                        (void*) &(*values.begin()),
                        values.size(),
                        SiloTraits<T>::datatype(),
                        optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBAddRegion
//------------------------------------------------------------------------------
inline
int
DBAddRegion(DBmrgtree_wrapper& tree,
            std::string reg_name,
            int info_bits,
            int max_children,
            std::string maps_name,
            std::vector<int>& seg_ids,
            std::vector<int>& seg_lens,
            std::vector<int>& seg_types,
            DBoptlist_wrapper& optlist) {
  int nsegs = seg_ids.size();
  VERIFY(seg_lens.size() == nsegs);
  VERIFY(seg_types.size() == nsegs);
  return DBAddRegion(tree.mDBmrgtree,
                     reg_name.c_str(),
                     info_bits,
                     max_children,
                     (maps_name.size() > 0 ? maps_name.c_str() : NULL),
                     nsegs,
                     (nsegs > 0 ? &(*seg_ids.begin()) : NULL),
                     (nsegs > 0 ? &(*seg_lens.begin()) : NULL),
                     (nsegs > 0 ? &(*seg_types.begin()) : NULL),
                     optlist.mOptlistPtr);
}

//------------------------------------------------------------------------------
// DBSetCwr
//------------------------------------------------------------------------------
inline
int
DBSetCwr(DBmrgtree_wrapper& tree,
         std::string path) {
  return DBSetCwr(tree.mDBmrgtree,
                  path.c_str());
}

//------------------------------------------------------------------------------
// DBGetCwr
//------------------------------------------------------------------------------
inline
const char*
DBGetCwr(DBmrgtree_wrapper& tree) {
  return DBGetCwr(tree.mDBmrgtree);
}

//------------------------------------------------------------------------------
// DBPutMrgtree
//------------------------------------------------------------------------------
inline
int
DBPutMrgtree(DBfile& file,
             std::string name,
             std::string mesh_name,
             DBmrgtree_wrapper& tree,
             DBoptlist_wrapper& optlist) {
  return DBPutMrgtree(&file,
                      name.c_str(),
                      mesh_name.c_str(),
                      tree.mDBmrgtree,
                      optlist.mOptlistPtr);
}

}

//------------------------------------------------------------------------------
// Typedef for vectors of optlists.
//------------------------------------------------------------------------------
typedef std::vector<silo::DBoptlist_wrapper*> vector_of_DBoptlist;

#endif
