//---------------------------------Spheral++----------------------------------//
// Slide Surface -- 
//----------------------------------------------------------------------------//

#include "FSISPH/SlideSurface.hh"
#include "DataBase/DataBase.hh"

#include <limits.h>
#include <vector>
#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------
template<typename Dimension>
SlideSurface<Dimension>::
SlideSurface(DataBase<Dimension>& dataBase,
             const std::vector<int> contactTypes):
  mIsActive(false),
  mNumNodeLists(0.0),
  mIsSlideSurface(){

    mNumNodeLists = dataBase.numNodeLists();

    // for our custom "map" (nodelisti,nodelistj) -> bool isSlide
    for(std::vector<int>::const_iterator it = contactTypes.begin();
        it != contactTypes.end();
        ++it){
      if (*it == 1){
        mIsActive=true;
        mIsSlideSurface.push_back(true);
      }else{
        mIsSlideSurface.push_back(false);
      }    
    }

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SlideSurface<Dimension>::
~SlideSurface() {
}

} // spheral namespace