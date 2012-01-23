#include <boost/python.hpp>
#include "Utilities/newtonRaphson.hh"

using namespace boost::python;

namespace Spheral {
//------------------------------------------------------------------------------
// Provide an overridable python interface to the Newton Raphson root finder.
//------------------------------------------------------------------------------
// A class we can override in python to provide the functor.
class NewtonRaphsonFunction {
public:
  NewtonRaphsonFunction() {};
  virtual ~NewtonRaphsonFunction() {};
  virtual std::pair<double, double> operator()(const double x) const { return std::pair<double, double>(x, 1.0); }
  virtual void thpt() {};
};

// The method we expose to python.
double newtonRaphsonFindRoot(const NewtonRaphsonFunction& functor,
                             float x1,
                             float x2,
                             const float xaccuracy = 1.0e-15,
                             const float yaccuracy = 1.0e-15,
                             const unsigned maxIterations = 100) {
  return Spheral::newtonRaphson(functor,
                                x1,
                                x2,
                                xaccuracy,
                                yaccuracy,
                                maxIterations);
}

//------------------------------------------------------------------------------
// Overridable wrapper for the functor class.
//------------------------------------------------------------------------------
struct NewtonRaphsonFunction_Wrapper: NewtonRaphsonFunction
{
    NewtonRaphsonFunction_Wrapper(PyObject* py_self_):
        NewtonRaphsonFunction(), py_self(py_self_) {}

    void thpt() {
        call_method< void >(py_self, "thpt");
    }

    void default_thpt() {
        NewtonRaphsonFunction::thpt();
    }

    PyObject* py_self;
};

//------------------------------------------------------------------------------
// Wrap the newtonRaphson method.
//------------------------------------------------------------------------------
void wrapNewtonRaphson() {
  class_<NewtonRaphsonFunction, 
         boost::noncopyable,
         NewtonRaphsonFunction_Wrapper>("NewtonRaphsonFunction", init<>())
    .def("thpt", &NewtonRaphsonFunction::thpt, &Spheral::NewtonRaphsonFunction_Wrapper::default_thpt)
    .def("__call__", &NewtonRaphsonFunction::operator())
    ;
  
  def("newtonRaphsonFindRoot", &newtonRaphsonFindRoot);
}

}
