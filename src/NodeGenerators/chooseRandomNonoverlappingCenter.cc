#include <random>

#include "chooseRandomNonoverlappingCenter.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Look for a position for the given shape inside a boundary and avoiding the
// given set of shapes.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
chooseRandomNonoverlappingCenter(typename Dimension::Vector& result,
                                 const typename Dimension::FacetedVolume& shape,
                                 const typename Dimension::FacetedVolume& boundary,
                                 const std::vector<typename Dimension::FacetedVolume>& existingShapes,
                                 const unsigned maxTries) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // A few useful constants.
  const Scalar r0 = Dimension::rootnu(shape.volume()/M_PI);
  const Scalar length = (boundary.xmax() - boundary.xmin()).maxElement();
  VERIFY(length > 0.0);
  const Vector xmin = boundary.xmin() - r0*Vector(1,1,1);

  // Prepare a random number generator.
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);

  // Can we find a place inside the boundary?
  unsigned iter = 0;
  bool done = false;
  while (iter < maxTries and not done) {
    ++iter;
    
    // Find a trial position in the allowed volume.  Note we gradually allow the center to be outside
    // the surface but touching as the iterations increase.
    result = xmin + (length + 2.0*r0)*Vector(dis(gen), dis(gen), dis(gen));
    while (not (boundary.contains(result) or
                (boundary.distance(result) < r0*double(iter)/maxTries))) {
      result = xmin + (length + 2.0*r0)*Vector(dis(gen), dis(gen), dis(gen));
    }
      
    // Check if the trial shape at this position would overlap any others.
    const FacetedVolume trialShape = shape + result;
    auto shapeItr = existingShapes.begin();
    while (shapeItr < existingShapes.end() and
           (not (trialShape.intersect(*shapeItr)))) ++shapeItr;
    done = (shapeItr == existingShapes.end());
  }

  // That's it.
  return iter;
}

}
