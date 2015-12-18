//------------------------------------------------------------------------------
// This is the C++ implementation of C interface functions for calling into
// Spheral from C/Fortran.
//------------------------------------------------------------------------------

extern "C" {
#include "SpheralC.h"
}

#include "Distributed/Communicator.hh"
#include "SpheralPseudoScript2D.hh"
#include "SpheralPseudoScript3D.hh"

//------------------------------------------------------------------------------
// spheral_set_communicator
//------------------------------------------------------------------------------
#ifdef USE_MPI
void spheral_set_communicator(MPI_Comm* comm) {
  Spheral::Communicator::communicator(*comm);
}
#endif

//------------------------------------------------------------------------------
// spheral_initialize3d
//------------------------------------------------------------------------------
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
                          const double xmax_z) {
  Spheral::SpheralPseudoScript3D::initialize((ASPH == 1),
                                             (XSPH == 1),
                                             (compatibleEnergy == 1),
                                             (vGradCorrection == 1),
                                             (hGradCorrection == 1),
                                             (sumMassDensity == 1),
                                             (useVelocityDt == 1),
                                             (ScalarQ == 1),
                                             (addDistributedBoundary == 1),
                                             (useDamage == 1),
                                             kernelType,
                                             piKernelType,
                                             gradKernelType,
                                             nmats,
                                             CFL,
                                             hmin,
                                             hmax,
                                             hminmaxratio,
                                             nPerh,
                                             Clinear,
                                             Cquadratic,
                                             xmin_x,
                                             xmin_y,
                                             xmin_z,
                                             xmax_x,
                                             xmax_y,
                                             xmax_z);
}

//------------------------------------------------------------------------------
// spheral_initialize2d
//------------------------------------------------------------------------------
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
                          const double xmax_y) {
  Spheral::SpheralPseudoScript2D::initialize((ASPH == 1),
                                             (XSPH == 1),
                                             (compatibleEnergy == 1),
                                             (vGradCorrection == 1),
                                             (hGradCorrection == 1),
                                             (sumMassDensity == 1),
                                             (useVelocityDt == 1),
                                             (ScalarQ == 1),
                                             (addDistributedBoundary == 1),
                                             (useDamage == 1),
                                             kernelType,
                                             piKernelType,
                                             gradKernelType,
                                             nmats,
                                             CFL,
                                             hmin,
                                             hmax,
                                             hminmaxratio,
                                             nPerh,
                                             Clinear,
                                             Cquadratic,
                                             xmin_x,
                                             xmin_y,
                                             xmax_x,
                                             xmax_y);
}

//------------------------------------------------------------------------------
// spheral_initialize_step3d
//------------------------------------------------------------------------------
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
                                 const int* particleType) {
  return Spheral::SpheralPseudoScript3D::initializeStep(nintpermat,
                                                        npermat,
                                                        mass,
                                                        massDensity,
                                                        position_x,
                                                        position_y,
                                                        position_z,
                                                        specificThermalEnergy,
                                                        velocity_x,
                                                        velocity_y,
                                                        velocity_z,
                                                        Hfield_xx,
                                                        Hfield_xy,
                                                        Hfield_xz,
                                                        Hfield_yy,
                                                        Hfield_yz,
                                                        Hfield_zz,
                                                        pressure,
                                                        deviatoricStress_xx,
                                                        deviatoricStress_xy,
                                                        deviatoricStress_xz,
                                                        deviatoricStress_yy,
                                                        deviatoricStress_yz,
                                                        deviatoricStress_zz,
                                                        soundSpeed,
                                                        bulkModulus,
                                                        shearModulus,
                                                        yieldStrength,
                                                        plasticStrain,
                                                        damage,
                                                        particleType);
}

//------------------------------------------------------------------------------
// spheral_initialize_step2d
//------------------------------------------------------------------------------
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
                                 const int* particleType) {
  return Spheral::SpheralPseudoScript2D::initializeStep(nintpermat,
                                                        npermat,
                                                        mass,
                                                        massDensity,
                                                        position_x,
                                                        position_y,
                                                        specificThermalEnergy,
                                                        velocity_x,
                                                        velocity_y,
                                                        Hfield_xx,
                                                        Hfield_xy,
                                                        Hfield_yy,
                                                        pressure,
                                                        deviatoricStress_xx,
                                                        deviatoricStress_xy,
                                                        deviatoricStress_yy,
                                                        soundSpeed,
                                                        bulkModulus,
                                                        shearModulus,
                                                        yieldStrength,
                                                        plasticStrain,
                                                        damage,
                                                        particleType);
}

//------------------------------------------------------------------------------
// spheral_update_state3d
//------------------------------------------------------------------------------
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
                            const int* particleType) {
  Spheral::SpheralPseudoScript3D::updateState(mass,
                                              massDensity,
                                              position_x,
                                              position_y,
                                              position_z,
                                              specificThermalEnergy,
                                              velocity_x,
                                              velocity_y,
                                              velocity_z,
                                              Hfield_xx,
                                              Hfield_xy,
                                              Hfield_xz,
                                              Hfield_yy,
                                              Hfield_yz,
                                              Hfield_zz,
                                              pressure,
                                              deviatoricStress_xx,
                                              deviatoricStress_xy,
                                              deviatoricStress_xz,
                                              deviatoricStress_yy,
                                              deviatoricStress_yz,
                                              deviatoricStress_zz,
                                              soundSpeed,
                                              bulkModulus,
                                              shearModulus,
                                              yieldStrength,
                                              plasticStrain,
                                              damage,
                                              particleType);
}

//------------------------------------------------------------------------------
// spheral_update_state2d
//------------------------------------------------------------------------------
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
                            const int* particleType) {
  Spheral::SpheralPseudoScript2D::updateState(mass,
                                              massDensity,
                                              position_x,
                                              position_y,
                                              specificThermalEnergy,
                                              velocity_x,
                                              velocity_y,
                                              Hfield_xx,
                                              Hfield_xy,
                                              Hfield_yy,
                                              pressure,
                                              deviatoricStress_xx,
                                              deviatoricStress_xy,
                                              deviatoricStress_yy,
                                              soundSpeed,
                                              bulkModulus,
                                              shearModulus,
                                              yieldStrength,
                                              plasticStrain,
                                              damage,
                                              particleType);
}

//------------------------------------------------------------------------------
// spheral_evaluate_derivatives3d
//------------------------------------------------------------------------------
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
                                    double* qwork) {
  Spheral::SpheralPseudoScript3D::evaluateDerivatives(massDensitySum,
                                                      DmassDensityDt,
                                                      DvelocityDt_x,
                                                      DvelocityDt_y,
                                                      DvelocityDt_z,
                                                      DxDt,
                                                      DyDt,
                                                      DzDt,
                                                      DspecificThermalEnergyDt,
                                                      DvelocityDx_xx,
                                                      DvelocityDx_xy,
                                                      DvelocityDx_xz,
                                                      DvelocityDx_yx,
                                                      DvelocityDx_yy,
                                                      DvelocityDx_yz,
                                                      DvelocityDx_zx,
                                                      DvelocityDx_zy,
                                                      DvelocityDx_zz,
                                                      DHfieldDt_xx,
                                                      DHfieldDt_xy,
                                                      DHfieldDt_xz,
                                                      DHfieldDt_yy,
                                                      DHfieldDt_yz,
                                                      DHfieldDt_zz,
                                                      HfieldIdeal_xx,
                                                      HfieldIdeal_xy,
                                                      HfieldIdeal_xz,
                                                      HfieldIdeal_yy,
                                                      HfieldIdeal_yz,
                                                      HfieldIdeal_zz,
                                                      DdeviatoricStressDt_xx,
                                                      DdeviatoricStressDt_xy,
                                                      DdeviatoricStressDt_xz,
                                                      DdeviatoricStressDt_yy,
                                                      DdeviatoricStressDt_yz,
                                                      DdeviatoricStressDt_zz,
                                                      qpressure,
                                                      qwork);
}

//------------------------------------------------------------------------------
// spheral_evaluate_derivatives2d
//------------------------------------------------------------------------------
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
                                    double* qwork) {
  Spheral::SpheralPseudoScript2D::evaluateDerivatives(massDensitySum,
                                                      DmassDensityDt,
                                                      DvelocityDt_x,
                                                      DvelocityDt_y,
                                                      DxDt,
                                                      DyDt,
                                                      DspecificThermalEnergyDt,
                                                      DvelocityDx_xx,
                                                      DvelocityDx_xy,
                                                      DvelocityDx_yx,
                                                      DvelocityDx_yy,
                                                      DHfieldDt_xx,
                                                      DHfieldDt_xy,
                                                      DHfieldDt_yy,
                                                      HfieldIdeal_xx,
                                                      HfieldIdeal_xy,
                                                      HfieldIdeal_yy,
                                                      DdeviatoricStressDt_xx,
                                                      DdeviatoricStressDt_xy,
                                                      DdeviatoricStressDt_yy,
                                                      qpressure,
                                                      qwork);
}

//------------------------------------------------------------------------------
// spheral_add_boundary
//------------------------------------------------------------------------------
void spheral_add_boundary(char * bcname, int ndims,
			  double p0, double p1, double p2,
			  double n0, double n1, double n2) {
  if(ndims==2) {
     Spheral::SpheralPseudoScript2D::addBoundary(bcname, p0, p1, n0, n1);
  }
  else {
     Spheral::SpheralPseudoScript3D::addBoundary(bcname, p0, p1, p2, n0, n1, n2);
  }
}

//------------------------------------------------------------------------------
// spheral_iterate_Hfield3d
//------------------------------------------------------------------------------
void spheral_iterate_Hfield3d(double* Hfield_xx,
                              double* Hfield_xy,
                              double* Hfield_xz,
                              double* Hfield_yy,
                              double* Hfield_yz,
                              double* Hfield_zz,
                              const int maxIterations,
                              const double tolerance) {
  Spheral::SpheralPseudoScript3D::iterateHfield(Hfield_xx,
                                                Hfield_xy,
                                                Hfield_xz,
                                                Hfield_yy,
                                                Hfield_yz,
                                                Hfield_zz,
                                                maxIterations,
                                                tolerance);
}

//------------------------------------------------------------------------------
// spheral_iterate_Hfield2d
//------------------------------------------------------------------------------
void spheral_iterate_Hfield2d(double* Hfield_xx,
                              double* Hfield_xy,
                              double* Hfield_yy,
                              const int maxIterations,
                              const double tolerance) {
  Spheral::SpheralPseudoScript2D::iterateHfield(Hfield_xx,
                                                Hfield_xy,
                                                Hfield_yy,
                                                maxIterations,
                                                tolerance);
}

//------------------------------------------------------------------------------
// spheral_update_connectivity3d
//------------------------------------------------------------------------------
void spheral_update_connectivity3d() {
  Spheral::SpheralPseudoScript3D::updateConnectivity();
}

//------------------------------------------------------------------------------
// spheral_get_connectivity3d
//------------------------------------------------------------------------------
void spheral_get_connectivity3d(int*** numNeighbors,
                                int**** connectivity) {
  Spheral::SpheralPseudoScript3D::getConnectivity(numNeighbors, connectivity);
}

//------------------------------------------------------------------------------
// spheral_get_num_materials3d
//------------------------------------------------------------------------------
int spheral_get_num_materials3d() {
  return Spheral::SpheralPseudoScript3D::getNumMaterials();
}

//------------------------------------------------------------------------------
// spheral_get_num_nodes3d
//------------------------------------------------------------------------------
int* spheral_get_num_nodes3d() {
  return Spheral::SpheralPseudoScript3D::getNumNodes();
}

//------------------------------------------------------------------------------
// spheral_get_num_internal_nodes3d
//------------------------------------------------------------------------------
int* spheral_get_num_internal_nodes3d() {
  return Spheral::SpheralPseudoScript3D::getNumInternalNodes();
}

//------------------------------------------------------------------------------
// spheral_get_num_ghost_nodes3d
//------------------------------------------------------------------------------
int* spheral_get_num_ghost_nodes3d() {
  return Spheral::SpheralPseudoScript3D::getNumGhostNodes();
}

//------------------------------------------------------------------------------
// spheral_update_connectivity2d
//------------------------------------------------------------------------------
void spheral_update_connectivity2d() {
  Spheral::SpheralPseudoScript2D::updateConnectivity();
}

//------------------------------------------------------------------------------
// spheral_get_connectivity2d
//------------------------------------------------------------------------------
void spheral_get_connectivity2d(int*** numNeighbors,
                                int**** connectivity) {
  Spheral::SpheralPseudoScript2D::getConnectivity(numNeighbors, connectivity);
}

//------------------------------------------------------------------------------
// spheral_get_num_materials2d
//------------------------------------------------------------------------------
int spheral_get_num_materials2d() {
  return Spheral::SpheralPseudoScript2D::getNumMaterials();
}

//------------------------------------------------------------------------------
// spheral_get_num_nodes2d
//------------------------------------------------------------------------------
int* spheral_get_num_nodes2d() {
  return Spheral::SpheralPseudoScript2D::getNumNodes();
}

//------------------------------------------------------------------------------
// spheral_get_num_internal_nodes2d
//------------------------------------------------------------------------------
int* spheral_get_num_internal_nodes2d() {
  return Spheral::SpheralPseudoScript2D::getNumInternalNodes();
}

//------------------------------------------------------------------------------
// spheral_get_num_ghost_nodes2d
//------------------------------------------------------------------------------
int* spheral_get_num_ghost_nodes2d() {
  return Spheral::SpheralPseudoScript2D::getNumGhostNodes();
}

