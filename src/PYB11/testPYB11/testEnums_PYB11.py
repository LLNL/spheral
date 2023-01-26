"""
This is a test module for enum binding in PYB11.
"""

from PYB11Generator import *

# List the files we want to include.
PYB11includes = ['<iostream>']

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

    enum class Furniture { chair, bed, couch };
  };

  enum class Color { black, white, red, blue, yellow };
}

namespace Bspace {
  template<typename T1, typename T2>
  class B {
  public:
    B() { std::cerr << "B<T1, T2>()" << std::endl; }
    enum class Rodent { mouse, squirrel, gerbil };
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

    Furniture = PYB11enum(("chair", "bed", "couch"))

#-------------------------------------------------------------------------------
# B
#-------------------------------------------------------------------------------
@PYB11namespace("Bspace")
@PYB11template("T1", "T2")
class B:

    def pyinit(self):
        "Default constructor B<%(T1)s, %(T2)s>."

    Rodent = PYB11enum(("mouse", "squirrel", "gerbil"))

# B<int, double>
Bintdouble = PYB11TemplateClass(B, template_parameters=("int", "double"))

#-------------------------------------------------------------------------------
# Color
#-------------------------------------------------------------------------------
Color = PYB11enum(("black", "white", "red", "blue", "yellow"),
                  namespace="Aspace")

#-------------------------------------------------------------------------------
# Attributes.
#-------------------------------------------------------------------------------
the_answer = PYB11attr("42")
what = PYB11attr('py::cast("The world")')
