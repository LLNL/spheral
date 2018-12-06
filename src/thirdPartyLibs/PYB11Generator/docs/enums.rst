.. _enums:

=====
Enums
=====

C++ enums are handled in a fairly straightforward manner as discussed in the `pybind11 docs <https://pybind11.readthedocs.io/en/stable/classes.html#enumerations-and-internal-types>`_.  Suppose we have the following enum in C++:

.. code-block:: cpp

   // A collection of adorable rodents
   enum Rodents {
     mouse = 0,
     rat = 1,
     squirrel = 2,
     hamster = 3,
     gerbil = 4,
     capybara = 5
   };

PYB11 uses the special method ``PYB11enum`` to declare enums, directly corresponding to the pybind11 construct ``py::enum_``.  We can bind our enumeration of Rodents using::

  Rodents = PYB11enum(("mouse", "rat", "squirrel", "hamster", "gerbil", "capybara"),
                      doc="A collection of adorable rodents")

See :func:`PYB11enum` for the full set of options to ``PYB11enum``.

It is also straightforward to declare an enum type that is inside a class scope; if we have a C++ class with an enum like the following:

.. code-block:: cpp

  class Homestararmy {
    public:

    enum members {
      HomestarRunner = 0,
      StrongSad = 1,
      Homsar = 2,
      PaintingOfGuyWithBigKnife = 3,
      FrankBennedetto = 4,
    };
  };

We can bind this enum using PYB11Generator with::

  class Homestararmy:

      members = PYB11enum(("HomestarRunner", "StrongSad", "Homsar", "PaintingOfGuyWithBigKnife", "FrankBennedetto"))

