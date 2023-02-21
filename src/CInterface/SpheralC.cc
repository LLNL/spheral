//------------------------------------------------------------------------------
// This is the C++ implementation of C interface functions for calling into
// Spheral from C/Fortran.
//------------------------------------------------------------------------------

#include "SpheralC.h"

#include <vector>

#include "Geometry/Dimension.hh"
#include "Distributed/Communicator.hh"
#include "SpheralPseudoScript.hh"

//------------------------------------------------------------------------------
// spheral_set_communicator
//------------------------------------------------------------------------------
#ifdef USE_MPI
void spheral_set_communicator(MPI_Comm* comm) {
  Spheral::Communicator::communicator(*comm);
}
#endif

//------------------------------------------------------------------------------
// spheral_initialize
//------------------------------------------------------------------------------
void spheral_initialize(const int      ndims,
                        const int      axisym,
                        const int      CRK,
                        const int      ASPH,
                        const int      XSPH,
                        const int      compatibleEnergy,
                        const int      totalEnergy,
                        const int      vGradCorrection,
                        const int      hGradCorrection,
                        const int      densityUpdate,
                        const int      sumMassDensity,
                        const int      useVelocityDt,
                        const int      Qoption,
                        const int      distributedBoundary,
                        const int      kernelType,
                        const int      nbspline,
                        const int      numinterp,
                        const int      rkorder,
                        const int      rkvolume,
                        const int      damage,
                        const unsigned nmats,
                        const double   CFL,
                        const double   hmin,
                        const double   hmax,
                        const double   hminmaxratio,
                        const double   nPerh,
                        const double   Clinear,
                        const double   Cquadratic,
                        const double*  xmincoords,
                        const double*  xmaxcoords) {
  switch(ndims) {
  case 3:
    {
      typedef Spheral::Dim<3> Dimension;
      Dimension::Vector xmin(xmincoords[0], xmincoords[1], xmincoords[2]),
        xmax(xmaxcoords[0], xmaxcoords[1], xmaxcoords[2]);
      Spheral::SpheralPseudoScript<Dimension>::initialize(false,
                                                          CRK,
                                                          ASPH,
                                                          XSPH,
                                                          compatibleEnergy,
                                                          totalEnergy,
                                                          vGradCorrection,
                                                          hGradCorrection,
                                                          densityUpdate,
                                                          sumMassDensity,
                                                          useVelocityDt,
                                                          Qoption,
                                                          distributedBoundary,
                                                          kernelType,
                                                          nbspline,
                                                          numinterp,
                                                          rkorder,
                                                          rkvolume,
                                                          damage,
                                                          nmats,
                                                          CFL,
                                                          hmin,
                                                          hmax,
                                                          hminmaxratio,
                                                          nPerh,
                                                          Clinear,
                                                          Cquadratic,
                                                          xmin,
                                                          xmax);
    }
    break;

  case 2:
    {
      typedef Spheral::Dim<2> Dimension;
      Dimension::Vector xmin(xmincoords[0], xmincoords[1]),
        xmax(xmaxcoords[0], xmaxcoords[1]);
      Spheral::SpheralPseudoScript<Dimension>::initialize(axisym,
                                                          CRK,
                                                          ASPH,
                                                          XSPH,
                                                          compatibleEnergy,
                                                          totalEnergy,
                                                          vGradCorrection,
                                                          hGradCorrection,
                                                          densityUpdate,
                                                          sumMassDensity,
                                                          useVelocityDt,
                                                          Qoption,
                                                          distributedBoundary,
                                                          kernelType,
                                                          nbspline,
                                                          numinterp,
                                                          rkorder,
                                                          rkvolume,
                                                          damage,
                                                          nmats,
                                                          CFL,
                                                          hmin,
                                                          hmax,
                                                          hminmaxratio,
                                                          nPerh,
                                                          Clinear,
                                                          Cquadratic,
                                                          xmin,
                                                          xmax);
    }
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_update_state
//------------------------------------------------------------------------------
void spheral_update_state(const int       ndims,
                          const unsigned* nintpermat,
                          const unsigned* npermat,
                          const double*   mass,
                          const double*   massDensity,
                          const double**  position,
                          const double*   specificThermalEnergy,
                          const double**  velocity,
                          const double**  Hfield,
                          const double*   pressure,
                          const double**  deviatoricStress,
                          const double*   soundSpeed,
                          const double*   bulkModulus,
                          const double*   shearModulus,
                          const double*   yieldStrength,
                          const double*   plasticStrain,
                          const double*   scalarDamage,
                          const int*      particleType,
                          const int*      regionNumber) {
  switch(ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::updateState(nintpermat,
                                                               npermat,
                                                               mass,
                                                               massDensity,
                                                               position,
                                                               specificThermalEnergy,
                                                               velocity,
                                                               Hfield,
                                                               pressure,
                                                               deviatoricStress,
                                                               soundSpeed,
                                                               bulkModulus,
                                                               shearModulus,
                                                               yieldStrength,
                                                               plasticStrain,
                                                               scalarDamage,
                                                               particleType,
                                                               regionNumber);
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::updateState(nintpermat,
                                                               npermat,
                                                               mass,
                                                               massDensity,
                                                               position,
                                                               specificThermalEnergy,
                                                               velocity,
                                                               Hfield,
                                                               pressure,
                                                               deviatoricStress,
                                                               soundSpeed,
                                                               bulkModulus,
                                                               shearModulus,
                                                               yieldStrength,
                                                               plasticStrain,
                                                               scalarDamage,
                                                               particleType,
                                                               regionNumber);
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_initialize_boundaries_and_physics
//------------------------------------------------------------------------------
void spheral_initialize_boundaries_and_physics(const int ndims) {
  switch(ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::initializeBoundariesAndPhysics();
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::initializeBoundariesAndPhysics();
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_initialize_step
//------------------------------------------------------------------------------
double spheral_initialize_step(const int ndims) {
  switch (ndims) {
  case 3:
    return Spheral::SpheralPseudoScript<Spheral::Dim<3>>::initializeStep();
    break;

  case 2:
    return Spheral::SpheralPseudoScript<Spheral::Dim<2>>::initializeStep();
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_dt_node
//------------------------------------------------------------------------------
int spheral_dt_node(const int ndims) {
  switch (ndims) {
  case 3:
    return Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtNode();
    break;

  case 2:
    return Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtNode();
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_dt_reason
//------------------------------------------------------------------------------
enum SpheralDtConstraint spheral_dt_reason(const int ndims) {
  switch (ndims) {
  case 3:
    if (Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "sound speed" ||
        Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "longitudinal sound speed" ||
        Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "deviatoric stress effective sound speed") {
      return SPHERAL_DT_Courant;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "artificial viscosity") {
      return SPHERAL_DT_Q;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "velocity divergence") {
      return SPHERAL_DT_Hydro;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "pairwise velocity difference" ||
             Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "velocity magnitude") {
      return SPHERAL_DT_Velocity;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<3>>::dtReason() == "acceleration") {
      return SPHERAL_DT_Accel;
    }
    else {
      return SPHERAL_DT_NoConstraint;
    }
    break;

  case 2:
    if (Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "sound speed" ||
        Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "longitudinal sound speed" ||
        Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "deviatoric stress effective sound speed") {
      return SPHERAL_DT_Courant;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "artificial viscosity") {
      return SPHERAL_DT_Q;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "velocity divergence") {
      return SPHERAL_DT_Hydro;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "pairwise velocity difference" ||
             Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "velocity magnitude") {
      return SPHERAL_DT_Velocity;
    }
    else if (Spheral::SpheralPseudoScript<Spheral::Dim<2>>::dtReason() == "acceleration") {
      return SPHERAL_DT_Accel;
    }
    else {
      return SPHERAL_DT_NoConstraint;
    }
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_evaluate_derivatives
//------------------------------------------------------------------------------
void spheral_evaluate_derivatives(const int ndims,
                                  double*   massDensitySum,
                                  double*   DmassDensityDt,
                                  double**  DvelocityDt,
                                  double**  DxDt,
                                  double*   DspecificThermalEnergyDt,
                                  double**  DvelocityDx,
                                  double**  DHfieldDt,
                                  double**  HfieldIdeal,
                                  double**  DdeviatoricStressDt,
                                  double*   qpressure,
                                  double*   qwork) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::evaluateDerivatives(massDensitySum,
                                                                       DmassDensityDt,
                                                                       DvelocityDt,
                                                                       DxDt,
                                                                       DspecificThermalEnergyDt,
                                                                       DvelocityDx,
                                                                       DHfieldDt,
                                                                       HfieldIdeal,
                                                                       DdeviatoricStressDt,
                                                                       qpressure,
                                                                       qwork);
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::evaluateDerivatives(massDensitySum,
                                                                       DmassDensityDt,
                                                                       DvelocityDt,
                                                                       DxDt,
                                                                       DspecificThermalEnergyDt,
                                                                       DvelocityDx,
                                                                       DHfieldDt,
                                                                       HfieldIdeal,
                                                                       DdeviatoricStressDt,
                                                                       qpressure,
                                                                       qwork);
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_add_boundary
//------------------------------------------------------------------------------
void spheral_add_boundary(const int     ndims,
                          const double* pcoords,
                          const double* ncoords) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3> >::addBoundary(Spheral::Dim<3>::Vector(pcoords[0], pcoords[1], pcoords[2]), 
                                                                Spheral::Dim<3>::Vector(ncoords[0], ncoords[1], ncoords[2]));
    break;
    
  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2> >::addBoundary(Spheral::Dim<2>::Vector(pcoords[0], pcoords[1]), 
                                                                Spheral::Dim<2>::Vector(ncoords[0], ncoords[1]));
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_periodic_boundary
//------------------------------------------------------------------------------
void spheral_periodic_boundary(const int     ndims,
                               const double* pcoords1,
                               const double* ncoords1,
                               const double* pcoords2,
                               const double* ncoords2) {
  switch(ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3> >::addPeriodicBoundary(Spheral::Dim<3>::Vector(pcoords1[0], pcoords1[1], pcoords1[2]), 
                                                                        Spheral::Dim<3>::Vector(ncoords1[0], ncoords1[1], ncoords1[2]),
                                                                        Spheral::Dim<3>::Vector(pcoords2[0], pcoords2[1], pcoords2[2]),
                                                                        Spheral::Dim<3>::Vector(ncoords2[0], ncoords2[1], ncoords2[2]));
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2> >::addPeriodicBoundary(Spheral::Dim<2>::Vector(pcoords1[0], pcoords1[1]), 
                                                                        Spheral::Dim<2>::Vector(ncoords1[0], ncoords1[1]),
                                                                        Spheral::Dim<2>::Vector(pcoords2[0], pcoords2[1]),
                                                                        Spheral::Dim<2>::Vector(ncoords2[0], ncoords2[1]));
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_iterate_Hfield
//------------------------------------------------------------------------------
void spheral_iterate_Hfield(const int    ndims,
                            double**     Hfield,
                            const int    maxIterations,
                            const double tolerance) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::iterateHfield(Hfield,
                                                                 maxIterations,
                                                                 tolerance);
    break;
    
  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::iterateHfield(Hfield,
                                                                 maxIterations,
                                                                 tolerance);
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_compute_fragments
//------------------------------------------------------------------------------
void spheral_compute_fragments(const int ndims,
                               double*   damage,
                               double    frag_radius,
                               double    frag_density,
                               double    frag_damage,
                               int*      fragments) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::computeFragmentID(damage,
                                                                     frag_radius,
                                                                     frag_density,
                                                                     frag_damage,
                                                                     fragments);
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::computeFragmentID(damage,
                                                                     frag_radius,
                                                                     frag_density,
                                                                     frag_damage,
                                                                     fragments);

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_sample_mesh
//------------------------------------------------------------------------------
void spheral_sample_mesh(const int      ndims,
                         const double*  xmincoords,
                         const double*  xmaxcoords,
                         const int*     nsamples,
                         double*        latticeDensity,
                         double**       latticeVelocity) {
  switch (ndims) {
  case 3:
    {
      typedef Spheral::Dim<3> Dimension;
      Dimension::Vector xmin(xmincoords[0], xmincoords[1], xmincoords[2]),
        xmax(xmaxcoords[0], xmaxcoords[1], xmaxcoords[2]);
      Spheral::SpheralPseudoScript<Dimension>::sampleLatticeMesh(xmin, xmax, nsamples,
                                                                 latticeDensity,
                                                                 latticeVelocity);
    }
    break;

  case 2:
    {    
      typedef Spheral::Dim<2> Dimension;
      Dimension::Vector xmin(xmincoords[0], xmincoords[1]), 
        xmax(xmaxcoords[0], xmaxcoords[1]);
      Spheral::SpheralPseudoScript<Dimension>::sampleLatticeMesh(xmin, xmax, nsamples,
                                                                 latticeDensity,
                                                                 latticeVelocity);
    }
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_polyhedral_mesh
//------------------------------------------------------------------------------
void spheral_polyhedral_mesh(const int      ndims,
                             int*           nnodes,
                             int*           nfaces,
                             int*           ncells,
                             double**       coords,
                             int**          facetonodes,
                             int**          nodecounts,
                             int**          celltofaces,
                             int**          facecounts,
                             int**          faceflags,
                             double**       volumes) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::polyhedralMesh(nnodes,
                                                                  nfaces,
                                                                  ncells,
                                                                  coords,
                                                                  facetonodes,
                                                                  nodecounts,
                                                                  celltofaces,
                                                                  facecounts,
                                                                  faceflags,
                                                                  volumes);
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::polyhedralMesh(nnodes,
                                                                  nfaces,
                                                                  ncells,
                                                                  coords,
                                                                  facetonodes,
                                                                  nodecounts,
                                                                  celltofaces,
                                                                  facecounts,
                                                                  faceflags,
                                                                  volumes);
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}


//------------------------------------------------------------------------------
// spheral_fill_volume
//------------------------------------------------------------------------------
void spheral_fill_volume(const int      ndims,
                         const int*     nnodes,
                         const int*     nfaces,
                         const double** coords,
                         const int*     conn,
                         const double   spacing,
                         const int      domain,
                         const int      ndomains,
                         double*        volume,
                         int*           nparticles,
                         double**       sphcoords) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::fillVolume(nnodes,
                                                              nfaces,
                                                              coords,
                                                              conn,
                                                              spacing,
                                                              domain,
                                                              ndomains,
                                                              volume,
                                                              nparticles,
                                                              sphcoords);
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::fillVolume(nnodes,
                                                              nfaces,
                                                              coords,
                                                              conn,
                                                              spacing,
                                                              domain,
                                                              ndomains,
                                                              volume,
                                                              nparticles,
                                                              sphcoords);
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_generate_cyl
//------------------------------------------------------------------------------
void spheral_generate_cyl(const int      ndims,
                          const int*     nnodes,
                          const double** coords,
                          const double** htensor,
                          const double** volume,
                          const double   frac,
                          int*           nparticles,
                          double**       sphcoords,
                          double**       sphhtensor,
                          double**       sphvolume) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::generateCylFromRZ(nnodes,
                                                                     coords,
                                                                     htensor,
                                                                     volume,
                                                                     frac,
                                                                     nparticles,
                                                                     sphcoords,
                                                                     sphhtensor,
                                                                     sphvolume);
    break;
    
  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::generateCylFromRZ(nnodes,
                                                                     coords,
                                                                     htensor,
                                                                     volume,
                                                                     frac,
                                                                     nparticles,
                                                                     sphcoords,
                                                                     sphhtensor,
                                                                     sphvolume);
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_update_connectivity
//------------------------------------------------------------------------------
void spheral_update_connectivity(const int ndims) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::updateConnectivity();
    break;

  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::updateConnectivity();
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_connectivity
//------------------------------------------------------------------------------
void spheral_get_connectivity(const int ndims,
                              int***    numNeighbors,
                              int****   connectivity) {
  switch (ndims) {
  case 3:
    Spheral::SpheralPseudoScript<Spheral::Dim<3>>::getConnectivity(numNeighbors, connectivity);
    break;
    
  case 2:
    Spheral::SpheralPseudoScript<Spheral::Dim<2>>::getConnectivity(numNeighbors, connectivity);
    break;

  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_materials
//------------------------------------------------------------------------------
int spheral_get_num_materials(const int ndims) {
  switch (ndims) {
  case 3:
    return Spheral::SpheralPseudoScript<Spheral::Dim<3>>::getNumMaterials();
    break;
    
  case 2:
    return Spheral::SpheralPseudoScript<Spheral::Dim<2>>::getNumMaterials();
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_nodes
//------------------------------------------------------------------------------
int* spheral_get_num_nodes(const int ndims) {
  switch (ndims) {
  case 3:
    return Spheral::SpheralPseudoScript<Spheral::Dim<3>>::getNumNodes();
    break;
    
  case 2:
    return Spheral::SpheralPseudoScript<Spheral::Dim<2>>::getNumNodes();
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_internal_nodes
//------------------------------------------------------------------------------
int* spheral_get_num_internal_nodes(const int ndims) {
  switch (ndims) {
  case 3:
    return Spheral::SpheralPseudoScript<Spheral::Dim<3>>::getNumInternalNodes();
    break;
    
  case 2:
    return Spheral::SpheralPseudoScript<Spheral::Dim<2>>::getNumInternalNodes();
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_ghost_nodes
//------------------------------------------------------------------------------
int* spheral_get_num_ghost_nodes(const int ndims) {
  switch (ndims) {
  case 3:
    return Spheral::SpheralPseudoScript<Spheral::Dim<3>>::getNumGhostNodes();
    break;

  case 2:
    return Spheral::SpheralPseudoScript<Spheral::Dim<2>>::getNumGhostNodes();
    break;
    
  default:
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}
