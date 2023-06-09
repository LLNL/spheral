Version vYYYY-MM-DD -- Release date YYYY-MM-DD

This release contains ...

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
