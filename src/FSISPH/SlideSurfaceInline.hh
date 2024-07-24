namespace Spheral{

//------------------------------------------------------------------------------
// set/get bool list of interactions 
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
isSlideSurface(const std::vector<bool> x) {
  mIsSlideSurface = x;
}
template<typename Dimension>
inline
std::vector<bool>
SlideSurface<Dimension>::
isSlideSurface() const {
  return mIsSlideSurface;
}

//------------------------------------------------------------------------------
// set/get bool if any slides exist
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
isActive(const bool x) {
  mIsActive = x;
}
template<typename Dimension>
inline
bool
SlideSurface<Dimension>::
isActive() const {
  return mIsActive;
}


//------------------------------------------------------------------------------
// set/get number of node lists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SlideSurface<Dimension>::
numNodeLists(const int x) {
  mNumNodeLists = x;
}
template<typename Dimension>
inline
int
SlideSurface<Dimension>::
numNodeLists() const {
  return mNumNodeLists;
}


//------------------------------------------------------------------------------
// more intelligable access
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool 
SlideSurface<Dimension>::
isSlideSurface(const int nodeListi, 
               const int nodeListj) const {
    const auto oneDimIndex = mNumNodeLists * nodeListi + nodeListj;
    return mIsSlideSurface[oneDimIndex];
}

//------------------------------------------------------------------------------
// this is from our old implementation where we would just reduce the AV
// based on the velocity and interface normal directions.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar 
SlideSurface<Dimension>::
slideCorrection(const typename Dimension::Scalar  smoothnessi,
                const typename Dimension::Scalar  smoothnessj, 
                const typename Dimension::Vector& normali,
                const typename Dimension::Vector& normalj,
                const typename Dimension::Vector& velocityi,
                const typename Dimension::Vector& velocityj) const {

    const auto ssij = this->pairwiseInterfaceSmoothness(smoothnessi,smoothnessj);
    const auto nij = this->pairwiseInterfaceNormal(smoothnessi,smoothnessj,normali,normalj);
    const auto vijhat = (velocityi-velocityj).unitVector();
    const auto fij = std::abs(nij.dot(vijhat));
    const auto slideCorr = (1.0-ssij) + (ssij)*fij*fij; 

    return slideCorr;      

}

//------------------------------------------------------------------------------
// weighted version of the slide correction
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar 
SlideSurface<Dimension>::
weightedSlideCorrection(const typename Dimension::Scalar  smoothnessi,
                        const typename Dimension::Scalar  smoothnessj, 
                        const typename Dimension::Vector& normali,
                        const typename Dimension::Vector& normalj,
                        const typename Dimension::Vector& velocityi,
                        const typename Dimension::Vector& velocityj,
                        const typename Dimension::Scalar  weighti,
                        const typename Dimension::Scalar  weightj) const {

    const auto ssij = this->weightedPairwiseInterfaceSmoothness(smoothnessi,smoothnessj,weighti,weightj);
    const auto nij = this->weightedPairwiseInterfaceNormal(smoothnessi,smoothnessj,normali,normalj,weighti,weightj);
    const auto vijhat = (velocityi-velocityj).unitVector();
    const auto fij = std::abs(nij.dot(vijhat));
    const auto slideCorr = (1.0-ssij) + (ssij)*fij*fij; 

    return slideCorr;      

}

//------------------------------------------------------------------------------
// return 1 if pairwise interaction is fully slide and ramp down to zero
// based on the max and min smoothness values for the interaction. These numbers
// were pulled out of a hat.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar 
SlideSurface<Dimension>::
pairwiseInterfaceSmoothness(const typename Dimension::Scalar smoothnessi,
                            const typename Dimension::Scalar smoothnessj) const {

    const auto ssij = (smoothnessi + smoothnessj)*0.5;
    return ( 1.0 - 10.0*std::min(std::max(0.95-(ssij),0.0),0.10) );
    //const auto ssMax = std::max(smoothnessi,smoothnessj);
    //const auto ssMin = std::min(smoothnessi,smoothnessj);

    //const auto ssijMax = ( 1.0 -  10.0*std::min(std::max(0.97-ssMax,0.0),0.10) ); // ramps down 0.97->0.87
    //const auto ssijMin = ( 1.0 -   5.0*std::min(std::max(0.90-ssMin,0.0),0.20) );  // ramps down 0.90->0.70
    //return ssijMax*ssijMin;      

}

//------------------------------------------------------------------------------
// weighted pairwise smoothness 
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar 
SlideSurface<Dimension>::
weightedPairwiseInterfaceSmoothness(const typename Dimension::Scalar smoothnessi,
                                    const typename Dimension::Scalar smoothnessj,
                                    const typename Dimension::Scalar weighti,
                                    const typename Dimension::Scalar weightj) const {

    const auto tiny = std::numeric_limits<double>::epsilon();
    const auto ssij = (smoothnessi*weighti + smoothnessj*weightj)/std::max(weighti+weightj,tiny);
    return ( 1.0 - 10.0*std::min(std::max(0.95-(ssij),0.0),0.10) );      

}


//------------------------------------------------------------------------------
// returns pairwise surface normal
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector 
SlideSurface<Dimension>::
pairwiseInterfaceNormal(const typename Dimension::Scalar  smoothnessi,
                        const typename Dimension::Scalar  smoothnessj,
                        const typename Dimension::Vector& normali,
                        const typename Dimension::Vector& normalj) const {
    return this->weightedPairwiseInterfaceNormal(smoothnessi,smoothnessj,normali,normalj,1.0,1.0);      
}

//------------------------------------------------------------------------------
// weighted pairwise surface normal 
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector 
SlideSurface<Dimension>::
weightedPairwiseInterfaceNormal(const typename Dimension::Scalar  smoothnessi,
                                const typename Dimension::Scalar  smoothnessj,
                                const typename Dimension::Vector& normali,
                                const typename Dimension::Vector& normalj,
                                const typename Dimension::Scalar  weighti,
                                const typename Dimension::Scalar  weightj) const {

    return (smoothnessj*weightj*normalj - smoothnessi*weighti*normali).unitVector();     

}

}
