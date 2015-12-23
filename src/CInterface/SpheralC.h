#ifndef __Spheral_C_interface__
#define __Spheral_C_interface__

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef __cplusplus
extern "C" {
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
  ----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  spheral_set_communicator

  Set the MPI communicator that Spheral will use.

  Arguments:
  comm (MPI_Communicator) : the MPI_COMM object to use.
  ----------------------------------------------------------------------------*/
#ifdef USE_MPI
void spheral_set_communicator(MPI_Comm* comm);
#endif

/*------------------------------------------------------------------------------
  spheral_initialize

  This method should be called once at the beginning of a calculation to set
  up the basic simulation parameters for Spheral.

  Returns: void
  Arguments:
    ASPH                  : flag selecting SPH or ASPH (0=false/1=true)
    ScalarQ               : flag selecting scalar or tensor viscosity (0=false/1=true)
    addDistribtedBoundary : flag whether we should use Spheral's internal method
                            of computing parallel ghost nodes (0=false/1=true)
    nmats                 : number of materials (NodeLists) to build
    hmin                  : minimum allowed smoothing scale
    hmax                  : maximum allowed smoothing scale
    hminmaxratio          : minimum allowed ratio of hmin/hmax for a given nodes
                            H tensor (1 => round, SPH)
    nPerh                 : target number of nodes per smoothing scale for
                            ideal H algorithm
    Clinear               : linear coefficient for artificial viscosity
    Cquadratic            : quadratic coefficient for artificial viscosity
    xmin_{x,y,z}          : minimum (x,y,z) coordinates in simulation volume for
                            neighbor searching
    xmax_{x,y,z}          : maximum (x,y,z) coordinates in simulation volume for
                            neighbor searching
  ----------------------------------------------------------------------------*/
void spheral_initialize3d(const int ASPH,
                          const int XSPH,
                          const int compatibleEnergy,
                          const int vGradCorrection,
                          const int hGradCorrection,
                          const int sumMassDensity,
                          const int useVelocityDt,
                          const int ScalarQ,
                          const int addDistributedBoundary,
                          const int useDamage,
                          const int kernelType,
                          const int piKernelType,
                          const int gradKernelType,
                          const unsigned nmats,
                          const double CFL,
                          const double hmin,
                          const double hmax,
                          const double hminmaxratio,
                          const double nPerh,
                          const double Clinear,
                          const double Cquadratic,
                          const double xmin_x,
                          const double xmin_y,
                          const double xmin_z,
                          const double xmax_x,
                          const double xmax_y,
                          const double xmax_z);

void spheral_initialize2d(const int ASPH,
                          const int XSPH,
                          const int compatibleEnergy,
                          const int vGradCorrection,
                          const int hGradCorrection,
                          const int sumMassDensity,
                          const int useVelocityDt,
                          const int ScalarQ,
                          const int addDistributedBoundary,
                          const int useDamage,
                          const int kernelType,
                          const int piKernelType,
                          const int gradKernelType,
                          const unsigned nmats,
                          const double CFL,
                          const double hmin,
                          const double hmax,
                          const double hminmaxratio,
                          const double nPerh,
                          const double Clinear,
                          const double Cquadratic,
                          const double xmin_x,
                          const double xmin_y,
                          const double xmax_x,
                          const double xmax_y);

/*------------------------------------------------------------------------------
  spheral_initialize_step

  This method should be called at the beginning of a time step to establish the 
  set of meshless point data that Spheral will be looking at.

  Returns: a vote for the timestep given the specified state data.
  Arguments:
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
  ----------------------------------------------------------------------------*/
double spheral_initialize_step3d(const unsigned* nintpermat,
                                 const unsigned* npermat,
                                 const double* mass,
                                 const double* massDensity,
                                 const double* position_x,
                                 const double* position_y,
                                 const double* position_z,
                                 const double* specificThermalEnergy,
                                 const double* velocity_x,
                                 const double* velocity_y,
                                 const double* velocity_z,
                                 const double* Hfield_xx,
                                 const double* Hfield_xy,
                                 const double* Hfield_xz,
                                 const double* Hfield_yy,
                                 const double* Hfield_yz,
                                 const double* Hfield_zz,
                                 const double* pressure,
                                 const double* deviatoricStress_xx,
                                 const double* deviatoricStress_xy,
                                 const double* deviatoricStress_xz,
                                 const double* deviatoricStress_yy,
                                 const double* deviatoricStress_yz,
                                 const double* deviatoricStress_zz,
                                 const double* soundSpeed,
                                 const double* bulkModulus,
                                 const double* shearModulus,
                                 const double* yieldStrength,
                                 const double* plasticStrain,
                                 const double* damage,
                                 const int* regionNumber,
                                 const int* particleType);

double spheral_initialize_step2d(const unsigned* nintpermat,
                                 const unsigned* npermat,
                                 const double* mass,
                                 const double* massDensity,
                                 const double* position_x,
                                 const double* position_y,
                                 const double* specificThermalEnergy,
                                 const double* velocity_x,
                                 const double* velocity_y,
                                 const double* Hfield_xx,
                                 const double* Hfield_xy,
                                 const double* Hfield_yy,
                                 const double* pressure,
                                 const double* deviatoricStress_xx,
                                 const double* deviatoricStress_xy,
                                 const double* deviatoricStress_yy,
                                 const double* soundSpeed,
                                 const double* bulkModulus,
                                 const double* shearModulus,
                                 const double* yieldStrength,
                                 const double* plasticStrain,
                                 const double* damage,
                                 const int* regionNumber,
                                 const int* particleType);

/*------------------------------------------------------------------------------
  spheral_update_state

  Updates the fluid state variables during a step if so desired.  This simply
  replaces the state originally specified in spheral_initialize_step, if for
  instance the host code is doing a multi-stage time integration step.

  Returns:  void
  Arguments:  same as the similarly named arguments in spheral_initialize_step
  ----------------------------------------------------------------------------*/
void spheral_update_state3d(const double* mass,
                            const double* massDensity,
                            const double* position_x,
                            const double* position_y,
                            const double* position_z,
                            const double* specificThermalEnergy,
                            const double* velocity_x,
                            const double* velocity_y,
                            const double* velocity_z,
                            const double* Hfield_xx,
                            const double* Hfield_xy,
                            const double* Hfield_xz,
                            const double* Hfield_yy,
                            const double* Hfield_yz,
                            const double* Hfield_zz,
                            const double* pressure,
                            const double* deviatoricStress_xx,
                            const double* deviatoricStress_xy,
                            const double* deviatoricStress_xz,
                            const double* deviatoricStress_yy,
                            const double* deviatoricStress_yz,
                            const double* deviatoricStress_zz,
                            const double* soundSpeed,
                            const double* bulkModulus,
                            const double* shearModulus,
                            const double* yieldStrength,
                            const double* plasticStrain,
                            const double* damage,
                            const int* regionNumber,
                            const int* particleType);

void spheral_update_state2d(const double* mass,
                            const double* massDensity,
                            const double* position_x,
                            const double* position_y,
                            const double* specificThermalEnergy,
                            const double* velocity_x,
                            const double* velocity_y,
                            const double* Hfield_xx,
                            const double* Hfield_xy,
                            const double* Hfield_yy,
                            const double* pressure,
                            const double* deviatoricStress_xx,
                            const double* deviatoricStress_xy,
                            const double* deviatoricStress_yy,
                            const double* soundSpeed,
                            const double* bulkModulus,
                            const double* shearModulus,
                            const double* yieldStrength,
                            const double* plasticStrain,
                            const double* damage,
                            const int* regionNumber,
                            const int* particleType);

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
void spheral_evaluate_derivatives3d(double* massDensitySum,
                                    double* DmassDensityDt,
                                    double* DvelocityDt_x,
                                    double* DvelocityDt_y,
                                    double* DvelocityDt_z,
                                    double* DxDt,
                                    double* DyDt,
                                    double* DzDt,
                                    double* DspecificThermalEnergyDt,
                                    double* DvelocityDx_xx,
                                    double* DvelocityDx_xy,
                                    double* DvelocityDx_xz,
                                    double* DvelocityDx_yx,
                                    double* DvelocityDx_yy,
                                    double* DvelocityDx_yz,
                                    double* DvelocityDx_zx,
                                    double* DvelocityDx_zy,
                                    double* DvelocityDx_zz,
                                    double* DHfieldDt_xx,
                                    double* DHfieldDt_xy,
                                    double* DHfieldDt_xz,
                                    double* DHfieldDt_yy,
                                    double* DHfieldDt_yz,
                                    double* DHfieldDt_zz,
                                    double* HfieldIdeal_xx,
                                    double* HfieldIdeal_xy,
                                    double* HfieldIdeal_xz,
                                    double* HfieldIdeal_yy,
                                    double* HfieldIdeal_yz,
                                    double* HfieldIdeal_zz,
                                    double* DdeviatoricStressDt_xx,
                                    double* DdeviatoricStressDt_xy,
                                    double* DdeviatoricStressDt_xz,
                                    double* DdeviatoricStressDt_yy,
                                    double* DdeviatoricStressDt_yz,
                                    double* DdeviatoricStressDt_zz,
                                    double* qpressure,
                                    double* qwork);

void spheral_evaluate_derivatives2d(double* massDensitySum,
                                    double* DmassDensityDt,
                                    double* DvelocityDt_x,
                                    double* DvelocityDt_y,
                                    double* DxDt,
                                    double* DyDt,
                                    double* DspecificThermalEnergyDt,
                                    double* DvelocityDx_xx,
                                    double* DvelocityDx_xy,
                                    double* DvelocityDx_yx,
                                    double* DvelocityDx_yy,
                                    double* DHfieldDt_xx,
                                    double* DHfieldDt_xy,
                                    double* DHfieldDt_yy,
                                    double* HfieldIdeal_xx,
                                    double* HfieldIdeal_xy,
                                    double* HfieldIdeal_yy,
                                    double* DdeviatoricStressDt_xx,
                                    double* DdeviatoricStressDt_xy,
                                    double* DdeviatoricStressDt_yy,
                                    double* qpressure,
                                    double* qwork);

/*------------------------------------------------------------------------------
  spheral_add_boundary

  Takes as input the dimension of the problem, the normal direction, and a
  point in the plane on which the boundary condition will be applied

  Returns:  void
  Arguments:  bcname      : type of boundary condition (symmetry, rigid, axial)
              ndims       : number of dimensions in the problem
              (p0,p1,p2)  : point in the plane
              (n0,n1,n2)  : normal direction at point
  ----------------------------------------------------------------------------*/
void spheral_add_boundary(char * bcname, int ndims,
                          double p0, double p1, double p2,
                          double n0, double n1, double n2) ;

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
void spheral_iterate_Hfield3d(double* Hfield_xx,
                              double* Hfield_xy,
                              double* Hfield_xz,
                              double* Hfield_yy,
                              double* Hfield_yz,
                              double* Hfield_zz,
                              const int maxIterations,
                              const double tolerance);

void spheral_iterate_Hfield2d(double* Hfield_xx,
                              double* Hfield_xy,
                              double* Hfield_yy,
                              const int maxIterations,
                              const double tolerance);

/*------------------------------------------------------------------------------
  spheral_update_connectivity

  Update the connectivity between nodes using Spheral's internal neighbor
  finding.

  Returns:  void
  Arguments:  None
  ----------------------------------------------------------------------------*/
void spheral_update_connectivity3d();

void spheral_update_connectivity2d();

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
void spheral_get_connectivity3d(int*** numNeighbors,
                                int**** connectivity);

void spheral_get_connectivity2d(int*** numNeighbors,
                                int**** connectivity);

/*------------------------------------------------------------------------------
  spheral_get_num_materials

  Get the current number of materials (Spheral NodeLists).

  Arguments:  None
  ----------------------------------------------------------------------------*/
int spheral_get_num_materials3d();

int spheral_get_num_materials2d();

/*------------------------------------------------------------------------------
  spheral_get_num_nodes

  Get the current total number of (A)SPH nodes per material.

  Arguments:  None
  ----------------------------------------------------------------------------*/
int* spheral_get_num_nodes3d();

int* spheral_get_num_nodes2d();

/*------------------------------------------------------------------------------
  spheral_get_num_internal_nodes

  Get the current total number of internal (A)SPH nodes per material.

  Arguments:  None
  ----------------------------------------------------------------------------*/
int* spheral_get_num_internal_nodes3d();

int* spheral_get_num_internal_nodes2d();

/*------------------------------------------------------------------------------
  spheral_get_num_ghost_nodes

  Get the current total number of ghost (A)SPH nodes per material.

  Arguments:  None
  ----------------------------------------------------------------------------*/
int* spheral_get_num_ghost_nodes3d();

int* spheral_get_num_ghost_nodes2d();

#ifdef __cplusplus
}
#endif

#endif
