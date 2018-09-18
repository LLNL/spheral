"""
This is a test module for class binding with inheritance in PYB11 (with templated classes).
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
  template<typename T>
  class A {
  private:
    T val;
  public:
    A(): val() { std::cerr << "A()" << std::endl; }
    virtual ~A() { std::cerr << "~A()" << std::endl; }
    virtual T do_something() const { std::cerr << "A::do_something" << std::endl; return 2*val; }
    virtual T do_something_else() const { std::cerr << "A::do_something_else" << std::endl; return 5*val; }
    virtual T yet_another_method() const { std::cerr << "A::yet_another_method" << std::endl; return 42*val; }
    T getval() const { return val; }
    void setval(T inval) { val = inval; }
  };
}

namespace Bspace {
  template<typename T>
  class B: public Aspace::A<T> {
  public:
    B(): Aspace::A<T>() { std::cerr << "B()" << std::endl; }
    virtual ~B() { std::cerr << "~B()" << std::endl; }
    virtual T do_something() const override { std::cerr << "B::do_something" << std::endl; return 100*this->getval();}
  };
}
"""

#-------------------------------------------------------------------------------
# A
#-------------------------------------------------------------------------------
@PYB11namespace("Aspace")
@PYB11template("T")
class A:

    def pyinit(self):
        "Default constructor."

    @PYB11virtual
    @PYB11const
    def do_something(self):
        "A virtual do_something method."
        return "%(T)s"

    @PYB11virtual
    @PYB11const
    def do_something_else(self):
        "A virtual do_something_else method."
        return "%(T)s"

    @PYB11virtual
    @PYB11const
    def yet_another_method(self):
        "A virtual yet_another_method."
        return "%(T)s"

#-------------------------------------------------------------------------------
# B
#-------------------------------------------------------------------------------
@PYB11namespace("Bspace")
@PYB11template("T")
class B(A):

    def pyinit(self):
        "Default constructor."

    @PYB11virtual
    @PYB11const
    def do_something(self):
        "B override of base do_something method."
        return "%(T)s"

#-------------------------------------------------------------------------------
# instantiations
#-------------------------------------------------------------------------------
Aint = PYB11TemplateClass(A, template_parameters = "int")
Bint = PYB11TemplateClass(B, template_parameters = "int")
