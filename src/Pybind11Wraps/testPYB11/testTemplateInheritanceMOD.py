"""
This is a test module for class binding with inheritance in PYB11 (with templated classes).
"""

from PYB11Generator import *

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
    virtual void yipes() = 0;
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
    virtual void yipes() { std::cerr << "B::yipes!" << std::endl; }
  };


  template<typename Ta, typename Tb>
  void some_function(const Ta& a, const Tb& b) { std::cerr << "some_function(" << a << " " << b << ")" << std::endl; }
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

    @PYB11pure_virtual
    def yipes(self):
        "Yipes!"
        return "void"

    @PYB11ignore
    @PYB11const
    def getval(self):
        return "%(T)s"

    @PYB11ignore
    def setval(self, inval = "%(T)s"):
        return "void"

    val = property(getval, setval, doc="Totally set the val.")

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

    @PYB11virtual
    def yipes(self):
        "Yipes!"
        return "void"

#-------------------------------------------------------------------------------
# some_function
#-------------------------------------------------------------------------------
@PYB11template("Ta", "Tb")
@PYB11namespace("Bspace")
def some_function(a = "const %(Ta)s&",
                  b = "const %(Tb)s&"):
    "A function of %(Ta)s and %(Tb)s."
    return "void"

#-------------------------------------------------------------------------------
# instantiations
#-------------------------------------------------------------------------------
Aint = PYB11TemplateClass(A, template_parameters = "int")
Bint = PYB11TemplateClass(B, template_parameters = "int")

some_function = PYB11TemplateFunction(some_function, template_parameters=("int", "double"))
