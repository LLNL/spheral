"""
This is a test module for class binding with inheritance in PYB11.
"""

from PYB11Generator import *

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
    virtual void do_something_else() const { std::cerr << "A::do_something_else" << std::endl; }
    virtual int yet_another_method() const { std::cerr << "A::yet_another_method" << std::endl; return 42; }
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
        return "void"

    @PYB11virtual
    @PYB11const
    def do_something_else(self):
        "A virtual do_something_else method."
        return "void"

    @PYB11virtual
    @PYB11const
    def yet_another_method(self):
        "A virtual yet_another_method."
        return "int"

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
        return "void"
