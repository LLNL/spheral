#-------------------------------------------------------------------------------
# RKCoefficients
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class RKCoefficients:
    "Carries the RK correction coefficients around"

    PYB11typedefs = """
  typedef RKCoefficients<%(Dimension)s> SelfType;
"""


    #...........................................................................
    # Comparators
    def __eq__(self):
        return

    def __ne__(self):
        return

    def __lt__(self):
        return

    def __gt__(self):
        return

    def __le__(self):
        return

    def __ge__(self):
        return

    #...........................................................................
    # Sequence methods
    @PYB11cppname("size")
    @PYB11const
    def __len__(self):
        return "unsigned"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference_internal")
    @PYB11implementation('[](SelfType& self, int i) { const int n = self.size(); if (i >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }, py::keep_alive<0,1>()')
    def __getitem__(self):
        return

    @PYB11implementation("[](SelfType& self, int i, const double v) { const int n = self.size(); if (i >= n) throw py::index_error(); self[(i %% n + n) %% n] = v; }")
    def __setitem__(self):
        "Set a value"

    @PYB11implementation("[](const SelfType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration"

    #...........................................................................
    # String representation
    @PYB11implementation("[](const SelfType& self) { std::ostringstream os; os << self; return os.str(); }")
    def __repr__(self):
        return

    #...........................................................................
    correctionOrder = PYB11readwrite(doc="The correction order of the coefficients")
    coeffs = PYB11readwrite(doc="The coefficients vector")
