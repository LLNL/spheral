Version vYYYY.MM.p -- Release date YYYY-MM-DD
==============================================
  * Important Notes:
    * External users of the code will need to supply config files for tpl-manager to find system libraries correctly. Steps to do this are detailed in the external user build guide. 

Notable changes include:

  * New features/ API changes:
    * Adding optional user specified smoothing scale method for SPH, FSISPH, and CRKSPH

  * Build changes / improvements:
    * tpl-manager.py will no longer use generic x86_64 configs for non LC systems. Users will be required to supply their own configs for pointing spack at external packages.

  * Bug Fixes / improvements:
    * Corrected an erroneous VERIFY in the P-alpha porosity constructor (with Fields of porosity and sound speed) that forced runs to stop even with correct input parameters
    * Fixed a bug in the standard ASPH hydros (ASPH, SolidASPH, and RZ varieties) that gave incorrect results.  FSI ad CRK models with ASPH smoothing scales were OK, but standard
      SPH using ASPH smoothing scales were simply incorrect for non-unit aspect ratio H's.  Also added ATS tests to help catch such errors going forward.

Version v2024.01.1 -- Release date 2024-02-17
==============================================
  * Important Notes:
    * This is a patch release for v2024.01.0.

Notable changes include:

  * New features/ API changes:
    * Adding an optional second-stage problem start-up hook to the Physics package interface: Physics::initializeProblemStartupDependencies.  The idea is to keep basic sizing
      of arrays and such in the first stage (Physics::initializeProblemStartup), while this new hook is used for updating any initial Physics state (and therefore provides a
      State and StateDerivatives object).

  * Build changes / improvements:
    * Improved the target export functionality.

  * Bug Fixes / improvements:
    * Fixed bug with ConstantBoundary in the presence of porosity with the new porosity models introduced in v2024.01.00.
    * Updating header lists for including Spheral modules in external projects.
    * Adding effective viscous pressure back to FSISPH.
    * Initial volumes for damage models were incorrectly not taking into account pore space when computing failure statistics for seeding flaws.  Fixed.

Version v2024.01.00 -- Release date 2024-01-19
==============================================
  * Important Notes:
    * The PolyClipper, BLT, and PYB11Generator submodules have been modified. Be sure to recursively update the submodules.  

Notable changes include:

  * New features/ API changes:
    * Adding P-alpha porosity model.
    * Updating treatment of various state variables in the presence of porosity.
    * Introduced a new common base class for porosity physics (PorosityModel), which PalphaPorosity and StrainPorosity share.
    * Revamped interaction UpdatePolicies with FieldLists:
      - UpdatePolicies have a new virtual method: clonePerField: True means when registering a FieldList copy the Policy for each Field in the FieldList; False means register the FieldList for update itself with the single instance of the Policy.
      - This change removes most of our redundant Field/FieldList update policies, and allows us to be more granular in applying different policies to single Field values in a FieldList.
    * Adding more Shadow Python interfaces wrapping our C++ classes, in particular PalphaPorosity and StrainPorosity.
    * EquationOfState now requires instances to provide \partial P/\partial \rho and \partial P/\partial \epsilon.  All current equations of state have been updated accordingly.
    * Tillotson and Gruneisen EOSs implementations updated a bit in the revamping.
    * Added more material options to MaterialPropertiesLib.py (mostly from Melosh's 89 book).

  * Build changes / improvements:
    * Spheral now provides First Class CMake support (using the BLT nomenclature). Spheral and its dependencies are now exported to simplify importing the project. To import Spheral into another project using CMake, use:
      ```
      find_package(Spheral_CXX <path_to_spheral_installation>)
      ```
    * CMake variables have a more consistent naming convention. Unused variables are removed.
    * Added ENABLE_DEV_BUILD option to improve build times during code development.
    * Upped our required C++ standard to 17.

  * Bug Fixes / improvements:
    * Fixed melt behavior in Steinberg-Guinan strength model, which was ignoring melt for damaged material.
    * Fixed range of dimensionless melt temperature for Johnson-Cook strength.
    * FSISPH new features and modifications to method. 
      * NOTE constructor inputs have changed.
      * strength implementation modified.
      * new features added including plane strain option and settable minP for interfaces.
      * new, more rigorous, interface and free surface tracking.
    * Fixed initialization of longitudinal sound speed and Youngs modulus for damage models.
    * Corrected some minor bugs/inconsistencies in the Tillotson EOS.
    * lcats updated to work with current TOSS4 machine configurations.
    * Updated various tests to make out automated testing more robust.

Version v2023-06-0 -- Release date 2023-06-20
==============================================

  * Important Notes:
    * This is the first release using Python 3.
    * Restart files are not backwards compatible with prior releases.

Notable changes include:

  * New features/ API changes:
    * We now use Python 3 rather than Python 2. The conversion from Python 2->3 is larger than most prior Python updates, which is part of why we held off for so long on this conversion. Python 2 is no longer supported however, and most other Python packages are only Python 3 now as well. Python comes with a tool (2to3) which can do most of the work converting a Python 2 script to Python 3, so we recommend when updating one of your old Python 2 run scripts or utilities you first let run 2to3 over that script.
    * FileIO has been simplified in how it handles Fields, which has reduced the API of FileIO.  However, this has two compromises:
      * Restart files are not backwards compatible with older releases.
      * Fields in restart files are pretty much binary blobs, and so are not easily inspected by hand.
    * We have removed support for the Python interface to Gnuplot, since the Python Gnuplot package had significant interface changes. Please convert such interactive graphics to Matplotlib, which has very similar support to our older Gnuplot utilities by importing SpheralMatplotlib.

  * Build changes / improvements:
    * PYB11Generator and PolyClipper are now submodules (in extern/)
    * Converted our PYB11 generator build rules to use the newly provided PYB11Generator Cmake rules.
    * toss_4_x86_64_ib system compatibility.
    * Updated spack version to v0.19.1
    * `--debug` and `--no-spec` options added to tpl-manager.py for outputing debug and skipping `spack spec` step.
    * `--init-only` option in tpl-manager.py will only initialize a local spack instance, skipping any TPL configuration.
    * TOSS4 compatibility for LC systems.
    * "risky" builds are installed on LC machines through gitlab CI.

  * Bug Fixes / improvements:
    * Fixed numerous compiler warnings with newer compilers such as G++ 9.4.
    * r-path for additional TPLs can be propogated to Spheral libraries with `SPHERAL_ADDITIONAL_RPATHS`.
    * The DEM package has received a significant updated.
      * Added simple analytic solid boundaries (planes,cylinder,sphere).
      * User given more control over DEM fast time stepping.
      * Simplified initialization and improved robustness.

Version v2023.03.0 -- Release date 2023-03-29
==============================================

This release contains ...

  * Important Notes:
    * Spheral now requires C++ 14
    * Spheral C++ libs are compiled into a single library (Spheral_CXX)

Notable changes include:

  * New features/ API changes:
    * New Discrete Element Model (DEM) physics package with linear-damped spring approach.
    * Adding a CUDA smoke test that can be called from the Spheral python API.
    * NVCC / CUDA 11 gitlab-ci jobs.
    * ATS default filters for non-MPI, debug and CUDA builds are injected into spheral-atstest script.
    * Latest Develop docker containers hosted on ghcr.io/llnl/spheral:latest.
    * External / offline builds are tested through github actions.
    * New polyhedral gravity solver.
    * Improved DEM timestep choice, so that points that do not overlap do not overly constrain the timestep.
    * Point potential gravity solver can now have non-unit metrics.
    * Pair-max damage algorithm now decouples points in different NodeLists or with different fragmentIDs.
    * Added direct support for FacetedVolumes to FileIO interface.
    * ANEOS now allows shallow copies to be made.
    * Polyhedron contain method can now optionally use Axom methods to test for containment.
    * Adding user specified functions for shear modulus and yield strength as a function of damage.

  * Build changes / improvements:
    * The C++ library interface is compiled into a single Spheral_CXX library. 
    * Previous Spheral C++ libraries are still CMAKE targets, but are "ojbect" libraries.
    * ATS bumped to version 7.0.9 for blueos smpi option support.
    * Eigen bumped to 3.4.0 for NVCC compatiblity.
    * C++ flag suppression is gaurded with build time CMake generators to only apply to C++ compilers.
    * Python runtime libraries are now managed through Spack / tpl-manager.
    * Added ENABLE_NAN_EXCEPTIONS (default OFF) Cmake flag to raise an exception when a NAN occurs (Gnu only).
    * Byte-compiling python installed in virtual spheral environment.
    * Invoking spheral no longer byte-compiles Python imported in a spheral script.
    * HDF5 links as a shared library. Static lib usage can be forced with ENABLE_STATIC_TPL.

  * Bug Fixes / improvements:
    * spheral-atstest scripts always point to locally installed ATS instance.
    * gitlab-ci report-results script for analyzing ATS CI runs.
    * Support for offline Spheral builds (provided TPLs are installed).
    * Fixes for restarting without regenerating the original node positions in the Python script.
    * Protected from division by zero in DEM when points coincide.
    * Corrected support for minimum pressure (intact and damaged) with porous materials.
    * Removed term driving damaged material to the reference density in solid hydros.
    * Added verbose flag to EquationOfState::specificThermalEnergyForPressure so users can see how the iterative search proceeds.

Version v2022.06.1 -- Release date 2022-06-24
=============================================

This is a bugfix release, which corrects a path problem that broke our convenient ANEOS
constructors using the provided input for quartz, dunite, and serpentine.


Version v2022.06.0 -- Release date 2022-06-09
=============================================

This release contains new features, bugfixes, and build improvements. Please see the
Spheral Documentation for more information about items in this release.

  * Important Notes: This Spheral release is the first release to support the new Spack based TPL (Third Party Library) management system. The previous CMake based TPL system has been deprecated and removed. The documentation has detailed instructions on how to use this new TPL system.

Notable changes include:

  * New features/ API changes:
    * Spack TPL management system.
    * Support for Gitlab-CI testing with Spheral.
    * Addition of SidreFileIO.
    * Support for Sidre Parallel IO (SPIO).
    * Full Documentation re-organization.
    * GSPH now has it's own pure-virtual GenericHydro class.
    * Adding some persistent state to the hydro objects to remeber criteria for diagnostics.
    * Moved inferace fields of SlideSurface class into the Hydro class.
    * FSISPH handles same-material damaged strength similarly to Spheral's default SPH solver.
    * New ProbabilisticDamageModel, which should be used in place of our prior Grady-Kipp implementations.
    * More Damage application options, and new defaults
    * Artificial viscosity for SPH variants now defaults to LimitedMonaghanGingoldViscosity.

  * Build changes / improvements:
    * Deleting CMake TPL system and all AutoTools Build system artifacts.
    * Removing support for mirrors in spack tpl system.

  * Bug Fixes / improvements:
    * Spheral fixed when running in Debug mode with MPI=Off.
    * Typos fixed in quickstart guide. https://github.com/LLNL/spheral/pull/116
    * Pedantic check for expired pointer to the RestartRegistrar. Ensures we don't call into deleted objects. 
    * Switching GammaLaw and Polytropic EOS to the isentrpic bulk modulus for consistency w/ Solid EOS.
    * Update scalar and tensor damage calc in FSISPH to be more consistent with SolidSPHHydro.
    * CullenDehnen segfault fix.

**Full Changelog**: https://github.com/LLNL/spheral/compare/2022.2.0-pre-spack...v2022.6.0

Version vYYYY.MM.p -- Release date YYYY-MM-DD
==============================================
  * Important Notes:

Notable changes include:

  * New features/ API changes:

  * Build changes / improvements:

  * Bug Fixes / improvements:
