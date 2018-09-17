"""
This is a test module for class binding with inheritance in PYB11.
"""

from PYB11Decorators import *
from PYB11STLmethods import *
from PYB11property import *
from PYB11class import *
from PYB11function import *

# List the files we want to include.
includes = ['<iostream>']

# We can specify arbitrary C++ to be inserted at the beginning of the file.
preamble = """
namespace Aspace {
  class A {
  public:
    A() { std::cerr << "A()" << std::endl; }
    virtual ~A() { std::cerr << "~A()" << std::endl; }
    virtual void do_something() const { std::cerr << "A::do_something" << std::endl; }
  };
}

namespace Bspace {
  class B: public Aspace::A {
  public:
    B(): Aspace::A() { std::cerr << "B()" << std::endl; }
    virtual ~B() { std::cerr << "~B()" << std::endl; }
    virtual void do_something() const override { std::cerr << "B::do_something" << std::endl; }
  };
}
"""

#-------------------------------------------------------------------------------
# A
#-------------------------------------------------------------------------------
@PYB11namespace("Aspace")
class A:

    def pyinit(self):
        "Default constructor."

    @PYB11virtual
    @PYB11const
    def do_something(self):
        "A virtual do_something method."
        return

#-------------------------------------------------------------------------------
# B
#-------------------------------------------------------------------------------
@PYB11namespace("Bspace")
class B(A):

    def pyinit(self):
        "Default constructor."

    @PYB11virtual
    @PYB11const
    def do_something(self):
        "B override of base do_something method."
        return
