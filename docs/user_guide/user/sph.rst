###################################################
What are these meshfree modeling methods?
###################################################

Spheral was conceived as a tool for the development and utilization of meshfree physics modeling algorithms, such as N-body gravitational models and Smoothed Particle Hydrodynamics (SPH).  In order to understand how best to use Spheral and it's applicability to a given problem, it is useful having a basic grasp of how these sorts of meshfree methods work.  This section is intended as a (very) brief introduction to these ideas -- the interested researcher is encouraged to dive into the detailed references to gain a deeper understanding.

SPH began in the astrophysics community which was looking for a method of modeling hydrodynamics coupled with N-body gravity methods.  Since N-body is particle based, it was natural to look for a hydrodynamics method which could also be particle based, and so SPH was born in the late 1970's.  SPH is indeed a method that solves the Lagrangian hydrodynamic conservation equations (mass, momentum, and energy1) on points that move freely about the problem.  Ordinary hydrodynamic methods ususally discretize space using a mesh of some sort, which carves space up into cells.  Meshfree methods instead distribute the mass of the problem amongst a finite number of points, which should be laid down conformally with the mass distribution being modeled.  As an example, consider a 2D cartoon example where we want represent a bounded circle of fluid using SPH:

.. subfigure::
   :layout-sm: AB
   :layout-lg: AB
   :layout-xl: AB
   :layout-xxl: AB
   :subcaptions: below
   :class-grid: outline

   .. image:: Circle.png
      :width: 98%
      :alt: Bounded circle of fluid.

   .. image:: Circle_SPH.png
      :width: 98%
      :alt: Evenly spaced SPH nodes to represent fluid

   Example of bounded circle of fluid in 2D we want to discretize using SPH

In this example we are representing the continuous green fluid on the left by evenly distributing a collection of SPH points in the volume of that fluid.  Each of the SPH points carries the physics variables we need for the physics problem we're trying to solve: mass, velocity, energy, temperature, deviatoric stress, etc.  SPH like methods take a discrete representation from a finite number of points such as this and extend that to a continuous representation by convolving a so-called interpolation kernel over the values on these points.  Interpolation kernels generally look like Gaussian functions of distance, except we generally use functions that have "compact support", which simply means that these functions fall to zero and terminate at some finite distance from their central point (unlike Gaussians which are formally infinite).  A cartoon representation of this interpolation process can be seen here:

.. figure:: SPH_sample_cartoon.*
   :width: 80%

   Notional SPH interpolation kernel centered on the red node.  Blue node have non-zero values for the kernel, while the black points do not and therefore do not contribute to the state of the blue point in question.

In this example the interpolation kernel is centered on the red point, and has non-zero values extending over the blue points.  However, it formally falls to zero outside this range, and therefore the black points do not contribute to the interpolation about the red point in question.

In SPH the interpolation kernel is represented by a function :math:`W(x^\alpha - x_i^\alpha, h)` where :math:`x^\alpha - x_i^\alpha` is the vector displacement between the sampling position and a point denoted by the index :math:`i`, and :math:`h` is the so-called smoothing scale which has units of length.  For most interpolation kernels this form can be reduced to a function of the nornalized distance :math:`W(\eta)`, where :math:`\eta^\alpha \equiv (x^\alpha - x_i^\alpha)/h` is the dimesionless coordinate vector and :math:`\eta = (\eta^\alpha \eta^\alpha)^{1/2}` is its magnitude.  (Note for mathematics throughout this guide we use the `Einstein summation convention <https://en.wikipedia.org/wiki/Einstein_notation>`_ for repeated indices.)  Common examples of functions that might be used as interpolation kernels include

  - Gaussian: :math:`W(\eta) = A \exp(-\eta^2)`
  - Wendland C4: :math:`W(\eta) = A \left(1 - \eta\right)^6 \left(1 + 6 \eta + \frac{35}{3} \eta^2\right), \forall \; \eta \le 1.0`

where the constant :math:`A` is used to enforce a volume normalization on the integral of :math:`W` such that :math:`\int W(\eta) \, dV = 1`.  Using this convention we can represent this volume convolution for SPH interpolation for a spatial field :math:`F(x^\alpha)` as

.. math::

   \langle F(x^\alpha) \rangle                &=       \int F(\prime{x}^\alpha) W(\prime{x}^\alpha - x^\alpha, h) dV \approx \sum_j V_j F(x_j^\alpha) W(x_j^\alpha - x^\alpha, h) \\
   \langle \partial_\beta F(x^\alpha) \rangle &\approx \sum_j V_j F(x_j^\alpha) \partial_\beta W(x_j^\alpha - x^\alpha, h) \\

So in the discrete approximation SPH in its simplest form provides numerical estimates of fields and their spatial gradients at any point in space (most crucially at the interpolation points themselves).  This same mathematical framework allows us to perform this spatial convolution over general partial differential equations (PDE's) and arrive at numerical approximations such as these for those PDE's as simple sums over the points near a given particle as functions of the interpolation kernel.

This leads us to a sometimes subtle but important distinction about these sorts of schemes: despite "Particle" being right there in the name of the method, SPH and its ilk are not really particle methods.  The points in SPH are best viewed as moving centers of interpolation, on which we are solving partial differential equations (PDE's), very similarly to how more traditional meshed methods such as finite-volume or finite-elements treat equations.  For this reason in Spheral we try to refer to use the term "nodes" rather than particles to refer to these points.  To be more concrete, in SPH we are solving the standard Lagragian conservation equations for mass, momentum, and energy either in the fluid or solid regimes:

Fluid equations:

.. math::

   \frac{D\rho}{Dt}        &= -\rho \partial_\alpha v^\alpha \\
   \frac{Dv^\alpha}{Dt}    &= -\rho^{-1} \partial_\alpha P \\
   \frac{D\varepsilon}{Dt} &= -\rho^{-1} P \partial_\alpha v^\alpha \\

Solid equations:

.. math::

   \frac{D\rho}{Dt}        &= -\rho \partial_\alpha v^\alpha \\
   \frac{Dv^\alpha}{Dt}    &= \rho^{-1} \partial_\beta \sigma^{\alpha \beta} \\
   \frac{D\varepsilon}{Dt} &= -\rho^{-1} \sigma^{\alpha \beta} \partial_\alpha v^\beta \\

where the standard fluid variables are

==========================================================================   =========================
:math:`\rho`                                                                 mass density             
:math:`V`                                                                    volume
:math:`v^\alpha`                                                             velocity vector          
:math:`P`                                                                    pressure                 
:math:`\varepsilon`                                                          specific thermal energy  
:math:`S^{\alpha \beta}`                                                     deviatoric stress        
:math:`\sigma^{\alpha \beta} = S^{\alpha \beta} - P \delta^{\alpha \beta}`   stress tensor
==========================================================================   =========================

These equations are solved in an SPH formalism by integrating the SPH interpolation 
