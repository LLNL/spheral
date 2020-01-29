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
                        const int      ScalarQ,
                        const int      distributedBoundary,
                        const int      kernelType,
                        const int      piKernelType,
                        const int      gradKernelType,
                        const int      nbspline,
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
  if (ndims == 3) {
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
                                                        ScalarQ,
                                                        distributedBoundary,
                                                        kernelType,
                                                        piKernelType,
                                                        gradKernelType,
                                                        nbspline,
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
  } else if (ndims == 2) {
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
                                                        ScalarQ,
                                                        distributedBoundary,
                                                        kernelType,
                                                        piKernelType,
                                                        gradKernelType,
                                                        nbspline,
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
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_initialize_step
//------------------------------------------------------------------------------
double spheral_initialize_step(const int       ndims,
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
                               const int*      particleType) {
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    return Spheral::SpheralPseudoScript<Dimension>::initializeStep(nintpermat,
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
                                                                   particleType);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    return Spheral::SpheralPseudoScript<Dimension>::initializeStep(nintpermat,
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
                                                                   particleType);
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_update_state
//------------------------------------------------------------------------------
void spheral_update_state(const int      ndims,
                          const double*  mass,
                          const double*  massDensity,
                          const double** position,
                          const double*  specificThermalEnergy,
                          const double** velocity,
                          const double** Hfield,
                          const double*  pressure,
                          const double** deviatoricStress,
                          const double*  soundSpeed,
                          const double*  bulkModulus,
                          const double*  shearModulus,
                          const double*  yieldStrength,
                          const double*  plasticStrain,
                          const double*  scalarDamage,
                          const int*     particleType) {
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::updateState(mass,
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
                                                         particleType);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::updateState(mass,
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
                                                         particleType);
  } else {
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
  if (ndims == 3) {
    // Build the vector packed types.
    typedef Spheral::Dim<3> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::evaluateDerivatives(massDensitySum,
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
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::evaluateDerivatives(massDensitySum,
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
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_add_boundary
//------------------------------------------------------------------------------
void spheral_add_boundary(const int     ndims,
                          const double* pcoords,
                          const double* ncoords) {
  if (ndims == 3) {
    Spheral::SpheralPseudoScript<Spheral::Dim<3> >::addBoundary(Spheral::Dim<3>::Vector(pcoords[0], pcoords[1], pcoords[2]), 
                                                                Spheral::Dim<3>::Vector(ncoords[0], ncoords[1], ncoords[2]));
  } else if (ndims == 2) {
    Spheral::SpheralPseudoScript<Spheral::Dim<2> >::addBoundary(Spheral::Dim<2>::Vector(pcoords[0], pcoords[1]), 
                                                                Spheral::Dim<2>::Vector(ncoords[0], ncoords[1]));
  } else {
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
  if (ndims == 3) {
    Spheral::SpheralPseudoScript<Spheral::Dim<3> >::addPeriodicBoundary(Spheral::Dim<3>::Vector(pcoords1[0], pcoords1[1], pcoords1[2]), 
                                                                        Spheral::Dim<3>::Vector(ncoords1[0], ncoords1[1], ncoords1[2]),
                                                                        Spheral::Dim<3>::Vector(pcoords2[0], pcoords2[1], pcoords2[2]),
                                                                        Spheral::Dim<3>::Vector(ncoords2[0], ncoords2[1], ncoords2[2]));
  } else if (ndims == 2) {
    Spheral::SpheralPseudoScript<Spheral::Dim<2> >::addPeriodicBoundary(Spheral::Dim<2>::Vector(pcoords1[0], pcoords1[1]), 
                                                                        Spheral::Dim<2>::Vector(ncoords1[0], ncoords1[1]),
                                                                        Spheral::Dim<2>::Vector(pcoords2[0], pcoords2[1]),
                                                                        Spheral::Dim<2>::Vector(ncoords2[0], ncoords2[1]));
  } else {
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
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::iterateHfield(Hfield,
                                                           maxIterations,
                                                           tolerance);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::iterateHfield(Hfield,
                                                           maxIterations,
                                                           tolerance);
  } else {
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
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::computeFragmentID(damage,
                                                               frag_radius,
                                                               frag_density,
                                                               frag_damage,
                                                               fragments);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::computeFragmentID(damage,
                                                               frag_radius,
                                                               frag_density,
                                                               frag_damage,
                                                               fragments);
  }
}

//------------------------------------------------------------------------------
// spheral_sample_mesh
//------------------------------------------------------------------------------
void spheral_sample_mesh(const int      ndims,
                         const double*  xmincoords,
                         const double*  xmaxcoords,
                         const int*     nsamples,
                         double*        latticeDensity) {
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    Dimension::Vector xmin(xmincoords[0], xmincoords[1], xmincoords[2]),
                      xmax(xmaxcoords[0], xmaxcoords[1], xmaxcoords[2]);
    Spheral::SpheralPseudoScript<Dimension>::sampleLatticeMesh(xmin, xmax, nsamples,
                                                               latticeDensity);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Dimension::Vector xmin(xmincoords[0], xmincoords[1]), 
                      xmax(xmaxcoords[0], xmaxcoords[1]);
    Spheral::SpheralPseudoScript<Dimension>::sampleLatticeMesh(xmin, xmax, nsamples,
                                                               latticeDensity);
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
                             int**          celltofaces) {
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::polyhedralMesh(nnodes,
                                                            nfaces,
                                                            ncells,
                                                            coords,
                                                            facetonodes,
                                                            celltofaces);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::polyhedralMesh(nnodes,
                                                            nfaces,
                                                            ncells,
                                                            coords,
                                                            facetonodes,
                                                            celltofaces);
  } else {
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
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::fillVolume(nnodes,
                                                        nfaces,
                                                        coords,
                                                        conn,
                                                        spacing,
                                                        domain,
                                                        ndomains,
                                                        volume,
                                                        nparticles,
                                                        sphcoords);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::fillVolume(nnodes,
                                                        nfaces,
                                                        coords,
                                                        conn,
                                                        spacing,
                                                        domain,
                                                        ndomains,
                                                        volume,
                                                        nparticles,
                                                        sphcoords);
  } else {
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
  if (ndims == 3) {
    typedef Spheral::Dim<3> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::generateCylFromRZ(nnodes,
                                                               coords,
                                                               htensor,
                                                               volume,
                                                               frac,
                                                               nparticles,
                                                               sphcoords,
                                                               sphhtensor,
                                                               sphvolume);
  } else if (ndims == 2) {
    typedef Spheral::Dim<2> Dimension;
    Spheral::SpheralPseudoScript<Dimension>::generateCylFromRZ(nnodes,
                                                               coords,
                                                               htensor,
                                                               volume,
                                                               frac,
                                                               nparticles,
                                                               sphcoords,
                                                               sphhtensor,
                                                               sphvolume);
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_update_connectivity
//------------------------------------------------------------------------------
void spheral_update_connectivity(const int ndims) {
  if (ndims == 3) {
    Spheral::SpheralPseudoScript<Spheral::Dim<3> >::updateConnectivity();
  } else if (ndims == 2) {
    Spheral::SpheralPseudoScript<Spheral::Dim<2> >::updateConnectivity();
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_connectivity
//------------------------------------------------------------------------------
void spheral_get_connectivity(const int ndims,
                              int***    numNeighbors,
                              int****   connectivity) {
  if (ndims == 3) {
    Spheral::SpheralPseudoScript<Spheral::Dim<3> >::getConnectivity(numNeighbors, connectivity);
  } else if (ndims == 2) {
    Spheral::SpheralPseudoScript<Spheral::Dim<2> >::getConnectivity(numNeighbors, connectivity);
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_materials
//------------------------------------------------------------------------------
int spheral_get_num_materials(const int ndims) {
  if (ndims == 3) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<3> >::getNumMaterials();
  } else if (ndims == 2) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<2> >::getNumMaterials();
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_nodes
//------------------------------------------------------------------------------
int* spheral_get_num_nodes(const int ndims) {
  if (ndims == 3) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<3> >::getNumNodes();
  } else if (ndims == 2) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<2> >::getNumNodes();
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_internal_nodes
//------------------------------------------------------------------------------
int* spheral_get_num_internal_nodes(const int ndims) {
  if (ndims == 3) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<3> >::getNumInternalNodes();
  } else if (ndims == 2) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<2> >::getNumInternalNodes();
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}

//------------------------------------------------------------------------------
// spheral_get_num_ghost_nodes
//------------------------------------------------------------------------------
int* spheral_get_num_ghost_nodes(const int ndims) {
  if (ndims == 3) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<3> >::getNumGhostNodes();
  } else if (ndims == 2) {
    return Spheral::SpheralPseudoScript<Spheral::Dim<2> >::getNumGhostNodes();
  } else {
    VERIFY2(false, "Error in SpheralC -- incorrect number of dimensions " << ndims << " requested.");
  }
}
