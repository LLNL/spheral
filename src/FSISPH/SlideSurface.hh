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

template<typename Dimension> class DataBase;

template<typename Dimension>
class SlideSurface {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  public:

    SlideSurface(DataBase<Dimension>& dataBase,
                 const std::vector<int> contactTypes);

    virtual ~SlideSurface();

    // return true if material iteraction is slide
    bool isSlideSurface(const int nodeListi, 
                        const int nodeListj) const;

    Scalar slideCorrection(const Scalar  smoothnessi,
                           const Scalar  smoothnessj, 
                           const Vector& normali,
                           const Vector& normalj,
                           const Vector& velocityi,
                           const Vector& velocityj) const;
    
    Scalar weightedSlideCorrection(const Scalar  smoothnessi,
                                   const Scalar  smoothnessj, 
                                   const Vector& normali,
                                   const Vector& normalj,
                                   const Vector& velocityi,
                                   const Vector& velocityj,
                                   const Scalar  weighti,
                                   const Scalar  weightj) const;

    Scalar pairwiseInterfaceSmoothness(const Scalar smoothnessi,
                                       const Scalar smoothnessj) const;

    Scalar weightedPairwiseInterfaceSmoothness(const Scalar smoothnessi,
                                               const Scalar smoothnessj,
                                               const Scalar weighti,
                                               const Scalar weightj) const;

    Vector pairwiseInterfaceNormal(const Scalar  smoothnessi,
                                   const Scalar  smoothnessj,
                                   const Vector& normali,
                                   const Vector& normalj) const;

    Vector weightedPairwiseInterfaceNormal(const Scalar  smoothnessi,
                                           const Scalar  smoothnessj,
                                           const Vector& normali,
                                           const Vector& normalj,
                                           const Scalar  weighti,
                                           const Scalar  weightj) const;

    
    std::vector<bool> isSlideSurface() const;
    void isSlideSurface(const std::vector<bool> x);

    bool isActive() const;
    void isActive(const bool x);

    int numNodeLists() const;
    void numNodeLists(const int x);


  private:
    //--------------------------- Private Interface ---------------------------//
    
    bool mIsActive;                                  // if there are 1 or more slide surface activate physics package
    int mNumNodeLists;                               // number of total node lists
    std::vector<bool> mIsSlideSurface;               // true if slide interaction between nodelists index --> numNodeList*nodeListi + nodeListj 

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
