//------------------------------------------------------------------------------
// GeometryRegistrar
//
// A singleton class for maintaining global information about the type of
// coordinate Spheral is running in.
//
// Created by JMO, Thu Mar 11 16:41:33 PST 2021
//------------------------------------------------------------------------------
#ifndef __Spheral_GeometryRegistrar__
#define __Spheral_GeometryRegistrar__

namespace Spheral {

class GeometryRegistrar {
public:
  //--------------------------- Public Interface ---------------------------//
  enum class CoordinateType {
    Cartesian = 0,
    Spherical = 1,
    RZ = 2,
  };

  // Get the instance.
  static GeometryRegistrar& instance();

  // The attribute we hang onto defining the local coordinate system.
  static CoordinateType coords;

private:
  //--------------------------- Private Interface --------------------------//
  static GeometryRegistrar* mInstancePtr;

  // No public constructors, destructor, or assignment.
  GeometryRegistrar();
  GeometryRegistrar(const GeometryRegistrar&);
  GeometryRegistrar& operator=(const GeometryRegistrar&);
  ~GeometryRegistrar();
};

}

#endif
