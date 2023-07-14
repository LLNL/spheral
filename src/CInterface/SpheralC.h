#ifndef __Spheral_C_interface__
#define __Spheral_C_interface__

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef SPHERAL_STATIC
#define SPHERALDLL_API
#else
#ifdef SPHERALDLL_EXPORTS
#define SPHERALDLL_API __declspec ( dllexport )
#else
#define SPHERALDLL_API __declspec ( dllimport )
#endif
#endif
#else
#define SPHERALDLL_API
#endif
/*------------------------------------------------------------------------------
  This is the main header for wrapping a C API calling into Spheral
  appropriate for use from C or Fortran.  This functionality is exposed here as
  a set of standalone C functions.

  Note, we should make every effort to ensure this interface remains
  pure C at this level!  Should make linking with C & Fortran simpler.

  A word about memory layout.  Spheral generally uses a NodeList per material,
  so the following interface allows the user to specify the number of materials
  which translates into that same number of NodeLists in Spheral.  If desired
  you can flatten all materials into a single material for passing into Spheral
  here, but you will lose the property that Spheral applies special rules for
  interactions between mateirals (in velocity gradients, mass density, strength
  interactions, etc).

  For arrays of node properties, we assume that properties are flattened across
  however many materials are specified into a single array, however.  So for
  the mass for instance, if we have 3 materials with (n1, n2, n3) nodes
  for each of these materials, the mass array would be a single C array of
  doubles layed out in memory like:
  [m(1,1), m(1,2), ..., m(1,n1), m(2,0), ..., m(2,n2), ..., m(3,0), ..., m(3,n3)]

  Using Spheral C involves a couple of different stages:

  1. Initialization at startup
  SpheralC calls should be performed once at the beginning of a simulation in
  the order:

    spheral_add_boundary(...)
      optionally called any number of times to establish boundaries

    spheral_initialize(...)
      sets various option run time options and switches
      - note, boundaries cannot be added after this call

    spheral_update_state(...)
      set the initial particle values 

    spheral_initialize_boundaries_and_physics(...)
      one-time initializations after intial particle data is provided

  2. Call order each physics steps
  After the above one-time program startup methods have been called, the
  time derivatives for the state variables can be queried at any step with the
  following sequence of calls:

    spheral_update_state(...)
      set the input particle state for this step

    spheral_initialize_step(...)
      performs physics and boundary preparation
      - also returns a vote for an appropriate time step

    spheral_evaluate_derivatives(...)
      computes the time derivatives per free particle
  ----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  spheral_set_communicator

  Set the MPI communicator that Spheral will use.

  Arguments:
  comm (MPI_Communicator) : the MPI_COMM object to use.
  ----------------------------------------------------------------------------*/
#ifdef USE_MPI
SPHERALDLL_API   void spheral_set_communicator(MPI_Comm* comm);
#endif

/*------------------------------------------------------------------------------
  spheral_initialize

  This method should be called once at the beginning of a calculation to set
  up the basic simulation parameters for Spheral.

  Returns: void
  Arguments:
    ndims                 : number of dimensions
    axisym                : flag for axisymmetric version
    CRK                   : flag for CRK (conservative reproducing kernel)
    ASPH                  : flag selecting SPH or ASPH (0=false/1=true)
    XSPH                  : flag selecting XSPH, i.e., move points with average fluid motion (0=false/1=true)
    compatibleEnergy      : flag for using the compatible energy method (0=false/1=true)
    totalEnergy           : flag for using the total energy derivative (0=false/1=true)
    vGradCorrection       : flag to select linear velocity gradient correction (0=false/1=true)
    hGradCorrection       : flag to select the gradh corection term (0=false/1=true)
    densityUpdate         : flag to select mass density definition (0=SumDensity, 1=RigorousSumDensity, 2=HybridSumDensity,
                                                                    3=IntegrateDensity, 4=VoronoiCellDensity, 5=SumVoronoiCellDensity,
                                                                    6=CorrectedSumDensity)
    sumMassDensity        : sum mass density over all NodeLists (0=false/1=true)
    useVelocityDt         : use the velocity magnitude to control time step (0=false/1=true)
    Qoption               : flag selecting the artificial viscosity (0=MonaghanGingold, 1=LimitedMonaghanGingold,
                                                                     2=TensorMonaghanGingold)
    distributedBoundary   : type of distributed boundary / neighbor method (0 or 1 = NestedGrid, 2 = Tree)
    kernelType            : select the generic interpolation kernel (0=BSpline, 1=Gaussian, 2=PiGaussian)
    nbspline              : order of kernel (if using B splines)
    numinterp             : number of kernel interpolation points
    rkorder               : order of CRK correction
    rkvolume              : RK volume definition (0=MassOverDensity, 1=SumVolume, 2=VoronoiVolume, 3=HullVolume, 4=HVolume)
    damage                : flag to feed back damage to Spheral
    nmats                 : number of materials (NodeLists) to build
    CFL                   : CFL timestep multiplier
    hmin                  : minimum allowed smoothing scale
    hmax                  : maximum allowed smoothing scale
    hminmaxratio          : minimum allowed ratio of hmin/hmax for a given nodes
                            H tensor (1 => round, SPH)
    nPerh                 : target number of nodes per smoothing scale for
                            ideal H algorithm
    Clinear               : linear coefficient for artificial viscosity
    Cquadratic            : quadratic coefficient for artificial viscosity
    xmincoords            : minimum (x,y,z) coordinates in simulation volume for
                            neighbor searching
    xmaxcoords            : maximum (x,y,z) coordinates in simulation volume for
                            neighbor searching
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
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
                          const double*  xmaxcoords);

/*------------------------------------------------------------------------------
  spheral_update_state

  Updates Spheral's versions of the fluid state variables

  Returns:  void
  Arguments:
    ndims                       : number of dimensions
    nintpermat                  : number of internal nodes per material (NodeList)
    npermat                     : total number nodes per material (NodeList)
    mass                        : mass per node
    massDensity                 : mass density per node
    position_{x,y,z}            : position coordinates for nodes
    specificThermalEnergy       : specific thermal energy per node
    velocity_{x,y,z}            : velocity components for nodes
    Hfield_{xx,xy,xz,           : Components of the H tensor for each node
               yy,yz,
                  zz}
    pressure                    : pressure per node
    deviatoricStress_{xx,xy,xz, : Components of the deviatoric stress per node.
                         yy,yz,
                            zz}
    soundSpeed                  : sound speed per node
    bulkModulus                 : bulkModulus per node
    shearModulus                : shearModulus per node
    yieldStrength               : yield strength per node
    plasticStrain               : plastic strain per node
    scalarDamage                : scalar damage per node
    particleType                : flag for active (1) or inactive (0) points
    regionNumber                : integer for the region ID
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
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
                            const int*      regionNumber);

/*------------------------------------------------------------------------------
  spheral_initialize_boundaries_and_physics

  This method should be called once at the beginning of a calculation in the
  sequence:
  spheral_initialize(...)                    - setup initial switches and objects
  spheral_update_state(...)                  - copy initial particle state fields
  spheral_initialize_boundaries_and_physics  - this method
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  void spheral_initialize_boundaries_and_physics(const int ndims);

/*------------------------------------------------------------------------------
  spheral_initialize_step

  This method should be called at the beginning of a time step to establish the 
  set of meshless point data that Spheral will be looking at, but after
  spheral_update_state.

  The order of methods to be called for each derivative evaluation should be
  spheral_update_state(...)
  spheral_initialize_step(...)
  spheral_evaluate_derivatives(...)

  Returns: a vote for the timestep given the specified state data.
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  double spheral_initialize_step(const int       ndims);

/*------------------------------------------------------------------------------
  spheral_dt_node, spheral_dt_reason

  These methods should be called after spheral_initialize_step to get the
  node ID which is controlling the time step, and the reason.

  The order of methods to be called for each derivative evaluation should be
  spheral_initialize_step(...)
  spheral_dt_node(...)
  spheral_dt_reason(...)

  Returns: information on the node controlling the time step
  ----------------------------------------------------------------------------*/
SPHERALDLL_API
  int spheral_dt_node(const int       ndims);

enum SpheralDtConstraint
{
  SPHERAL_DT_NoConstraint = 0,
  SPHERAL_DT_Courant,
  SPHERAL_DT_Q,
  SPHERAL_DT_Hydro,
  SPHERAL_DT_Velocity,
  SPHERAL_DT_Accel
};

SPHERALDLL_API
  enum SpheralDtConstraint spheral_dt_reason(const int       ndims);

/*------------------------------------------------------------------------------
  spheral_evaluate_derivatives

  Evaluates the current derivatives (time and space) based on the state last
  handed in.

  Returns:  void
  Arguments:
    massDensitySum      : The computed sum of the density per node
    DmassDensityDt      : The time derivative for the mass density per node
    DvelocityDt_{x,y,z} : The components of the acceleration computed per node
    DspecificThermalEnergyDt : The computed time derivative of the specific
                               thermal energy per node
    DvelocityDx_{xx,xy,xz, : The components of the spatial velocity gradient
                 yx,yy,yz,   per node
                 zx,zy,zz}
    DHfieldDt_{xx,xy,xz,   : The components of the time derivative of the H
                  yy,yz,     tensor per node
                     zz}
    HfieldIdeal_{xx,xy,xz, : Components of the measured best fit H tensor
                    yy,yz,   for the current point positions per node
                       zz}
    DdeviatoricStressDt_{xx,xy,xz, : Components of the time derivative of the 
                            yy,yz,   deviatoric stress per node
                               zz}
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
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
                                    double*   qwork);

/*------------------------------------------------------------------------------
  spheral_add_boundary

  Takes as input the dimension of the problem, the normal direction, and a
  point in the plane on which the boundary condition will be applied

  Returns:  void
  Arguments:  ndims       : number of dimensions in the problem
              (p0,p1,p2)  : point in the plane
              (n0,n1,n2)  : normal direction at point
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
void spheral_add_boundary(const int     ndims,
                          const double* pcoords,
                          const double* ncoords);

/*------------------------------------------------------------------------------
  spheral_periodic_boundary

  Takes as input the dimension of the problem, the normal direction, and a
  point in the plane on which the boundary condition will be applied

  Returns:  void
  Arguments:  ndims       : number of dimensions in the problem
              (p0,p1,p2)  : point in the plane
              (n0,n1,n2)  : normal direction at point
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
void spheral_periodic_boundary(const int     ndims,
                               const double* pcoords1,
                               const double* ncoords1,
                               const double* pcoords2,
                               const double* ncoords2);

/*------------------------------------------------------------------------------
  spheral_iterate_Hfield

  Takes as input a guess for the H tensor field, and optimizes it based on the
  chosen algorithm set in initialize (SPH or ASPH).

  Returns:  void
  Arguments:  Hfield_{xx,xy,xz, : Components of the per node H tensor field.
                         yy,yz,   Starts with the users first guess, and at the
                            zz}   end contains the final best guess.
              maxIterations     : maximum allowed iterations to try and optimize
                                  the H field
              tolerance         : maximum relative tolerance for changes in any
                                  nodes H tensor, defining stop condition for
                                  the iteration
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  void spheral_iterate_Hfield(const int    ndims,
                              double**     Hfield,
                              const int    maxIterations,
                              const double tolerance);

/*------------------------------------------------------------------------------
  spheral_compute_fragments

  Takes as input the scalar damage, and computes the fragment ID field

  Returns:  void
  Arguments:  damage        : pointer to the scalar damage for each particle
              frag_density  : dust tolerance for density
              frag_damage   : dust tolerance for damage
              fragments     : pointer to the fragment ID number
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  void spheral_compute_fragments(const int ndims,
                                 double*   damage,
                                 double    frag_radius,
                                 double    frag_density,
                                 double    frag_damage,
                                 int*      fragments);

/*------------------------------------------------------------------------------
  spheral_sample_mesh

  Takes as input a bounding box and sample points in each direction,
  and returns SPH field values at each point

  Returns:  void
  Arguments:
  ----------------------------------------------------------------------------*/
SPHERALDLL_API
void spheral_sample_mesh(const int      ndims,
                         const double*  xmincoords,
                         const double*  xmaxcoords,
                         const int*     nsamples,
                         double*        latticeDensity,
                         double**       latticeVelocity);

/*------------------------------------------------------------------------------
  spheral_polyhedral_mesh

  Takes as input a set of SPH particles and returns the polyhedral cells

  Returns:  void
  Arguments:
  ----------------------------------------------------------------------------*/
SPHERALDLL_API
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
                             double**       volumes);

/*------------------------------------------------------------------------------
  spheral_fill_volume

  Takes as input a set of coordinates for a hex or tet mesh and
  returns an even distribution of SPH particle coordinates

  Returns:  void
  Arguments:
  ----------------------------------------------------------------------------*/
SPHERALDLL_API
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
                         double**       sphcoords);

/*------------------------------------------------------------------------------
  spheral_generate_cyl

  Takes as input a set of 2D axisymmetric coordinates and revolves it to
  create a 3D cylindrical distribution of particle coordinates

  Returns:  void
  Arguments:
  ----------------------------------------------------------------------------*/
SPHERALDLL_API
void spheral_generate_cyl(const int      ndims,
                          const int*     nnodes,
                          const double** coords,
                          const double** htensor,
                          const double** volume,
                          const double   frac,
                          int*           nparticles,
                          double**       sphcoords,
                          double**       sphhtensor,
                          double**       sphvolume);

/*------------------------------------------------------------------------------
  spheral_update_connectivity

  Update the connectivity between nodes using Spheral's internal neighbor
  finding.

  Returns:  void
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  void spheral_update_connectivity(const int ndims);

/*------------------------------------------------------------------------------
  spheral_get_connectivity

  Get the current connectivity information.

  Returns:  void
  Arguments:
    numNeighbors :
      A three deep integer array indexed as
        numNeighbors[imat][i][jmat]
      where imat and jmat are material indices, and i is the node index we're
      interested in.  So numNeighbors[imat][i][jmat] is the number of neighbor 
      points in material (jmat) for the node (i) in material (imat).
    connectivity :
      A four deep integer array indexed as
        result[imat][i][jmat][j]
      where (imat, i) are the material number and index of node we are querying
      neighbor connectivity for, and result[imat][i][jmat][j] is the j'th
      neighbor in material (jmat) for the i'th node in material (imat).
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  void spheral_get_connectivity(const int ndims,
                                int***    numNeighbors,
                                int****   connectivity);

/*------------------------------------------------------------------------------
  spheral_get_num_materials

  Get the current number of materials (Spheral NodeLists).
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  int spheral_get_num_materials(const int ndims);

/*------------------------------------------------------------------------------
  spheral_get_num_nodes

  Get the current total number of (A)SPH nodes per material.
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  int* spheral_get_num_nodes(const int ndims);

/*------------------------------------------------------------------------------
  spheral_get_num_internal_nodes

  Get the current total number of internal (A)SPH nodes per material.
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  int* spheral_get_num_internal_nodes(const int ndims);

/*------------------------------------------------------------------------------
  spheral_get_num_ghost_nodes

  Get the current total number of ghost (A)SPH nodes per material.
  ----------------------------------------------------------------------------*/
SPHERALDLL_API 
  int* spheral_get_num_ghost_nodes(const int ndims);

#ifdef __cplusplus
}
#endif

#endif
