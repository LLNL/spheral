Version vYYYY.MM.pp -- Release date 20YY-MM-DD
==============================================

This release contains ...

  * Important Notes:

Notable changes include:

  * New features/ API changes:
    * New Discrete Element Model (DEM) physics package with linear-damped spring approach
    
  * Build changes / improvements:

  * Bug Fixes / improvements:


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

Version v2022.06.1 -- Release date 2022-06-24
=============================================

This is a bugfix release, which corrects a path problem that broke our convenient ANEOS
constructors using the provided input for quartz, dunite, and serpentine.
