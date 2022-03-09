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

    SlideSurface(const std::vector<int> contactTypes);

    SlideSurface(const std::vector<bool> contactTypes);

    virtual ~SlideSurface();

    // Slide specific methods
    //------------------------------------------------------
    // return true if material iteraction is slide
    bool isSlideSurface(const int nodeListi, 
                        const int nodeListj) const;

    Scalar slideCorrection(const int nodeListi,
                           const int i, 
                           const int nodeListj,
                           const int j,
                           const Vector vi,
                           const Vector vj) const;

    // Physics Package methods
    //------------------------------------------------------
    //virtual 
    //TimeStepType dt(const DataBase<Dimension>& dataBase, 
    //                const State<Dimension>& state,
    //                const StateDerivatives<Dimension>& derivs,
    //                const Scalar currentTime) const  override;

    // Tasks we do once on problem startup.
    virtual
    void initializeProblemStartup(DataBase<Dimension>& dataBase);

    // calls calc func for surface normals
    virtual
    void initialize(const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    const StateDerivatives<Dimension>& derivs,
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
                    const Scalar /*time*/,
                    const Scalar /*dt*/,
                    const TableKernel<Dimension>& W);
    
    // override w/ non-op
    //virtual
    //void evaluateDerivatives(const Scalar time,
    //                         const Scalar dt,
    //                         const DataBase<Dimension>& dataBase,
    //                         const State<Dimension>& state,
    //                               StateDerivatives<Dimension>& derivatives) const override;

    // register the surface normal field
    virtual
    void registerState(DataBase<Dimension>& dataBase,
                       State<Dimension>& state);

    // override w/ non-op
    //virtual
    //void registerDerivatives(DataBase<Dimension>& dataBase,
    //                         StateDerivatives<Dimension>& derivs) override;


    const FieldList<Dimension, Vector>& surfaceNormals() const;
    const FieldList<Dimension, Scalar>& surfaceFraction() const;
    const FieldList<Dimension, Scalar>& surfaceSmoothness() const;

    std::vector<bool> isSlideSurface() const;
    void isSlideSurface(const std::vector<bool> x);

    bool isActive() const;
    void isActive(const bool x);

    int numNodeLists() const;
    void numNodeLists(const int x);

    // non-ops for now
    //virtual std::string label() const override { return "SlideSurface" ; }
    //virtual void dumpState(FileIO& file, const std::string& pathName) const;
    //virtual void restoreState(const FileIO& file, const std::string& pathName);

  private:
    //--------------------------- Private Interface ---------------------------//
    bool mIsActive;                                  // is there are 1 or more slide surface activate physics package
    int mNumNodeLists;                               // number of total node lists
    std::vector<bool> mIsSlideSurface;               // true if slide interaction between nodelists index --> numNodeList*nodeListi + nodeListj 
    
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
