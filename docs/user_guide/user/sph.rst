###################################################
What are these meshfree modeling methods?
###################################################

Spheral was conceived as a tool for the development and utilization of meshfree physics modeling algorithms, such as N-body gravitational models and Smoothed Particle Hydrodynamics (SPH).  In order to understand how best to use Spheral and it's applicability to a given problem, it is useful having a basic grasp of how these sorts of meshfree methods work.  This section is intended as a (very) brief introduction to these ideas -- the interested researcher is encouraged to dive into the detailed references to gain a deeper understanding.

SPH began in the astrophysics community which was looking for a method of modeling hydrodynamics coupled with N-body gravity methods.  Since N-body is particle based, it was natural to look for a hydrodynamics method which could also be particle based, and so SPH was born in the late 1970's.  


.. |circle| image:: Circle.*

.. |SPH_circle| image:: Circle_SPH.*

.. table:: Discretizing a circle of material in SPH
   :align: center

   +-----------------------------------------------------+------------------------------------------------------+
   |                     |circle|                        |                    |SPH_circle|                      |
   | Bounded circle of fluid we wish to represent in 2D. |      Evenly spaced SPH nodes to represent fluid      |
   +-----------------------------------------------------+------------------------------------------------------+

Hurgle gurgle

.. figure:: SPH_sample_cartoon.*
   :width: 80%

   Notional SPH interpolation kernel centered on the red node.  Blue node have non-zero values for the kernel, while the black points do not and therefore do not contribute to the state of the blue point in question.
