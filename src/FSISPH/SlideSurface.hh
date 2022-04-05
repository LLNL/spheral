//---------------------------------Spheral++----------------------------------//
// Slide Surface -- 
//----------------------------------------------------------------------------//
#ifndef __Spheral_SlideSurface_hh__
#define __Spheral_SlideSurface_hh__

#include <string>
#include <vector>
#include <utility>

#include "Geometry/Dimension.hh"

namespace Spheral {

enum class SurfaceNormalMethod {
  SameMaterialSurfaceNormals = 0,
  DifferentMaterialSurfaceNormals = 1,
  AllMaterialSurfaceNormals = 2,
  MassWeightedSurfaceNormals = 3,
};

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class Boundary;
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
class FileIO;

template<typename Dimension>
class SlideSurface {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename std::vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  public:

    SlideSurface(DataBase<Dimension>& dataBase,
                 const std::vector<int> contactTypes,
                 const SurfaceNormalMethod surfaceNormalMethod,
                 const bool normalAreSmoothed=false,
                 const bool gradientsAreCorrected=true);

    virtual ~SlideSurface();

    // return true if material iteraction is slide
    bool isSlideSurface(const int nodeListi, 
                        const int nodeListj) const;

    // old method spits out multiplier to reduce av at interfaces
    Scalar slideCorrection(const int nodeListi,
                           const int i, 
                           const int nodeListj,
                           const int j,
                           const Vector vi,
                           const Vector vj) const;

    Scalar pairwiseSurfaceSmoothness(const int nodeListi,
                                     const int i, 
                                     const int nodeListj,
                                     const int j) const;

    Vector pairwiseSurfaceNormal(const int nodeListi,
                                 const int i, 
                                 const int nodeListj,
                                 const int j) const;
    
    Vector weightedPairwiseSurfaceNormal(const int nodeListi,
                                         const int i, 
                                         const int nodeListj,
                                         const int j,
                                         const Scalar weighti,
                                         const Scalar weightj) const;

    virtual
    void initializeProblemStartup(DataBase<Dimension>& dataBase);

    virtual
    void initialize(const DataBase<Dimension>& dataBase,
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs,
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
                    const Scalar /*time*/,
                    const Scalar /*dt*/,
                    const TableKernel<Dimension>& W);

    virtual
    void registerState(DataBase<Dimension>& dataBase,
                       State<Dimension>& state);

    void computeSurfaceSmoothness(const TableKernel<Dimension>& W,
                                  const DataBase<Dimension>& dataBase,
                                        State<Dimension>& state);

    void computeSurfaceNormals(const TableKernel<Dimension>& W,
                               const DataBase<Dimension>& dataBase,
                                     State<Dimension>& state);
    
    void smoothSurfaceNormals(const TableKernel<Dimension>& W,
                              const DataBase<Dimension>& dataBase,
                                    State<Dimension>& state);

    const FieldList<Dimension, Vector>& surfaceNormals() const;
    const FieldList<Dimension, Scalar>& surfaceFraction() const;
    const FieldList<Dimension, Scalar>& surfaceSmoothness() const;

    std::vector<bool> isSlideSurface() const;
    void isSlideSurface(const std::vector<bool> x);

    bool isActive() const;
    void isActive(const bool x);

    bool normalsAreSmoothed() const;
    void normalsAreSmoothed(const bool x);

    bool gradientsAreCorrected() const;
    void gradientsAreCorrected(const bool x);

    int numNodeLists() const;
    void numNodeLists(const int x);

    SurfaceNormalMethod surfaceNormalMethod() const;
    void surfaceNormalMethod(const SurfaceNormalMethod method);

  private:
    //--------------------------- Private Interface ---------------------------//
    
    bool mIsActive;                                  // if there are 1 or more slide surface activate physics package
    bool mNormalsAreSmoothed;                        // add additional operation to smooth interface normals
    bool mGradientsAreCorrected;                     // do we want to do the M-correction

    int mNumNodeLists;                               // number of total node lists
    std::vector<bool> mIsSlideSurface;               // true if slide interaction between nodelists index --> numNodeList*nodeListi + nodeListj 
    SurfaceNormalMethod mSurfaceNormalMethod;        // what subset of nodes do we base the normal calc off of

    FieldList<Dimension, Vector> mSurfaceNormals;    // surface normals between nodelists     
    FieldList<Dimension, Scalar> mSurfaceFraction;   // fraction of dissimilar neighbor volume     
    FieldList<Dimension, Scalar> mSurfaceSmoothness; // smoothness metric (0-1)    
  
    SlideSurface();
    SlideSurface(const SlideSurface&);
    SlideSurface& operator=(const SlideSurface&);

};

}


#include "SlideSurfaceInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SlideSurface;
}

#endif
